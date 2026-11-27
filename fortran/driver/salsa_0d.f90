! SALSA++
! Copyright (c) 2026 Wageningen University & Research
! SPDX-License-Identifier: Apache-2.0
!
! 0D driver for the Fortran reference. It uses upstream's own LES glue on a 5 x 5 x 3 grid (one interior point):
! the state is set up by setPrognosticVariables/setDiagnosticVariables, run_SALSA computes the SALSA tendencies,
! and SALSA_diagnostics moves evaporated cloud and rain back to aerosol. The sequence follows UCLALES-SALSA:
! - initialisation (init.f90): one run_SALSA call with initialize = .true., tendencies applied with clipping at zero;
! - every time step (step.f90): run_SALSA, tend_constrain2, forward Euler update, SALSA_diagnostics;
!   here SALSA_diagnostics runs first, before run_SALSA, as SALSA++ cleans the state at the start of `exec`.
! Settings are read from the file NAMELIST (&salsa for SALSA, &driver for this driver); forcing and initial
! conditions from the NetCDF input file.
program salsa_0d
    use grid, only: nxp, nyp, nzp, level, dtlt
    use classFieldArray, only: FieldArray
    use mo_field_state, only: Prog, Diag, SALSA_tracers_4d
    use mo_submctl, only: spec, nbins, ncld, nprc, in1a, fn1a, in2a, fn2a, in2b, fn2b, pi6, mair, rg, mwa, aerobins
    use mo_salsa_types, only: aero
    use mo_salsa_init, only: define_salsa, salsa_initialize
    use mo_salsa_sizedist, only: size_distribution
    use mo_progn_state, only: setPrognosticVariables, a_rp, a_naerop, a_maerop, a_gaerop
    use mo_diag_state, only: setDiagnosticVariables, a_temp, a_press, a_dn, a_rsl, a_rh
    use mo_aux_state, only: setGridSpacings, aetot
    use mo_salsa_driver, only: run_SALSA
    use constrain_SALSA, only: SALSA_diagnostics, tend_constrain2
    use util, only: getMassIndex
    use salsa_0d_io, only: Input_data, read_input, create_output, write_output, close_output
    use salsa_0d_version, only: git_commit
    implicit none

    ! Grid point of the 5 x 5 x 3 grid where SALSA runs.
    integer, parameter :: kp = 2, ip = 3, jp = 3

    ! &driver settings.
    character(len=256) :: case_name = '', input_file = 'input.nc', output_file = 'output_fortran.nc'
    real :: dt = 1.
    real :: w = 0.
    integer :: nsteps = 400
    namelist /driver/ case_name, input_file, output_file, dt, nsteps, w

    type(Input_data) :: d
    type(FieldArray) :: axes
    real, allocatable, target :: a_sclrp(:,:,:,:), a_sclrt(:,:,:,:)
    real, allocatable :: wp(:,:,:), nactd(:,:,:,:), vactd(:,:,:,:)
    character(len=8) :: outputlist(1) = ['none']
    integer :: nwet, nscl, memsize, istp, unit

    open(newunit=unit, file='NAMELIST', status='old', action='read')
    read(unit, nml=driver)
    close(unit)
    dtlt = dt

    ! Reads &salsa from the file NAMELIST.
    call define_salsa(level)
    call salsa_initialize(1000., 0.)

    ! Size of the scalar block as laid out by setPrognosticVariables for level 4, isgstyp <= 1, no piggybacking
    ! and no charge tracers: tp, rp, number and mass (dry compounds + water) per bin, and 5 gas tracers.
    nwet = spec%getNSpec(type="wet")
    nscl = 2 + (nbins + ncld + nprc)*(1 + nwet) + 5

    allocate(a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
    a_sclrp = 0.
    a_sclrt = 0.

    allocate(wp(nzp,nxp,nyp), nactd(nzp,nxp,nyp,ncld), vactd(nzp,nxp,nyp,nwet*ncld))
    wp = w

    Prog = FieldArray()
    Diag = FieldArray()
    axes = FieldArray()
    memsize = 0
    call setPrognosticVariables(a_sclrp, a_sclrt, Prog, outputlist, level, 1, .false., nzp, nxp, nyp, nscl)
    call setDiagnosticVariables(Diag, outputlist, memsize, level, 0, .false., nzp, nxp, nyp)
    call Prog%getByGroup("SALSA_4d", SALSA_tracers_4d)

    ! Aerosol bin lower limits, used to place evaporated rain back into the aerosol bins.
    call setGridSpacings(axes, level, nzp, nxp, nyp)
    aetot%d(1:nbins) = aerobins(1:nbins)

    call read_input(trim(input_file), d)
    if (size(d%time) < nsteps + 1) error stop "salsa_0d: the input file has fewer than nsteps + 1 forcing times"

    ! Initial state.
    call set_thermodynamics(d%temperature(1), d%pressure(1))
    a_rp%d(kp,ip,jp) = d%relative_humidity(1)*a_rsl%d(kp,ip,jp)
    a_rh%d = a_rp%d/a_rsl%d
    call set_initial_aerosol(d)
    call set_initial_gases(d)

    ! Initialisation call, as in UCLALES-SALSA init.f90 (initialize = .true. equilibrates the aerosol water).
    call run_SALSA(Diag, Prog, nzp, nxp, nyp, nwet, wp, nactd, vactd, dt, 0., 0, level, .false., .true.)
    ! Apply the tendencies, clipping the SALSA tracers at zero (SALSAInit). Index 1 is tp, 2 is rp.
    a_sclrp(kp:nzp,ip:nxp-2,jp:nyp-2,1:2) = a_sclrp(kp:nzp,ip:nxp-2,jp:nyp-2,1:2) &
                                          + dt*a_sclrt(kp:nzp,ip:nxp-2,jp:nyp-2,1:2)
    a_sclrp(kp:nzp,ip:nxp-2,jp:nyp-2,3:) = max(a_sclrp(kp:nzp,ip:nxp-2,jp:nyp-2,3:) &
                                          + dt*a_sclrt(kp:nzp,ip:nxp-2,jp:nyp-2,3:), 0.)
    a_sclrt = 0.
    call SALSA_diagnostics(.true., .false.)
    a_rh%d = a_rp%d/a_rsl%d

    call create_output(trim(output_file), trim(case_name), trim(input_file), read_text_file('NAMELIST'), git_commit)
    call write_output(0.)

    ! Time loop, as in UCLALES-SALSA step.f90. Step istp sees the forcing at its end time istp*dt
    ! (in the LES, the dynamics update the state before SALSA); run_SALSA gets the time at the start of the step.
    do istp = 1, nsteps
        if (abs(d%time(istp+1) - istp*dt) > 1.e-6*dt) error stop "salsa_0d: input time does not match istp*dt"

        call set_thermodynamics(d%temperature(istp+1), d%pressure(istp+1))
        a_rh%d = a_rp%d/a_rsl%d
        call SALSA_diagnostics(.false., .false.)

        call run_SALSA(Diag, Prog, nzp, nxp, nyp, nwet, wp, nactd, vactd, dt, (istp-1)*dt, istp, level, .false., .false.)
        call tend_constrain2()
        a_sclrp(kp:nzp-1,ip:nxp-2,jp:nyp-2,:) = a_sclrp(kp:nzp-1,ip:nxp-2,jp:nyp-2,:) &
                                              + dt*a_sclrt(kp:nzp-1,ip:nxp-2,jp:nyp-2,:)
        a_sclrt = 0.

        call write_output(istp*dt)
    end do

    call close_output()
    print '(a,i0,a)', 'salsa_0d: case ' // trim(case_name) // ', ', nsteps, ' steps, output in ' // trim(output_file)

contains

    ! Temperature and pressure from the forcing; air density and saturation mixing ratio as in the old driver.f90.
    subroutine set_thermodynamics(temperature, pressure)
        real, intent(in) :: temperature, pressure
        a_temp%d = temperature
        a_press%d = pressure
        a_dn%d = mair*pressure/(rg*temperature)
        a_rsl%d = mwa*saturation_vapour_pressure(temperature)/(mair*pressure)
    end subroutine set_thermodynamics


    ! Saturation vapour pressure over liquid water (Pa): 6th-order polynomial in T - 273.15 K (Lowe, 1977),
    ! as satvaph2o in SALSA-standalone mo_salsa_dynamics.f90 (commit 41da7aa, Apache 2.0).
    real function saturation_vapour_pressure(temperature)
        real, intent(in) :: temperature
        real, parameter :: a0 = 6.107799961, a1 = 4.436518521e-1, a2 = 1.428945805e-2, a3 = 2.650648471e-4, &
                           a4 = 3.031240396e-6, a5 = 2.034080948e-8, a6 = 6.136820929e-11
        real :: tc
        tc = temperature - 273.15
        saturation_vapour_pressure = 100.*(a0 + a1*tc + a2*tc**2 + a3*tc**3 + a4*tc**4 + a5*tc**5 + a6*tc**6)
    end function saturation_vapour_pressure


    ! Initial aerosol from lognormal modes, as UCLALES-SALSA does for a uniform profile (isdtyp = 0):
    ! number from size_distribution for regimes 1a+2a and 2b, mass = fraction * number * bin centre volume * density.
    ! Regime 1a holds only SO4 and OC, in the ratio of their regime-a fractions.
    subroutine set_initial_aerosol(d)
        type(Input_data), intent(in) :: d
        real :: core(nbins), nsect_a(1,1,nbins), nsect_b(1,1,nbins), ndist(nbins)
        real :: vf_a(nbins), vf_b(nbins), frac_oc_1a
        integer :: ndry, s, j, iso4, ioc
        logical :: found

        ndry = spec%getNSpec(type="dry")

        ! Volume fractions of the active compounds, by name.
        vf_a = 0.
        vf_b = 0.
        do j = 1, size(d%species_dry_name)
            found = .false.
            do s = 1, ndry
                if (trim(spec%names(s)) == trim(d%species_dry_name(j))) then
                    vf_a(s) = d%vol_frac_a(j)
                    vf_b(s) = d%vol_frac_b(j)
                    found = .true.
                end if
            end do
            if (.not. found .and. (d%vol_frac_a(j) /= 0. .or. d%vol_frac_b(j) /= 0.)) &
                error stop "salsa_0d: initial volume fraction given for a compound that is not in listspec"
        end do

        core = pi6*aero(1,1,1:nbins)%dmid**3
        nsect_a = 0.
        nsect_b = 0.
        call size_distribution(1, 1, 1, size(d%number_a), in1a, fn2a, d%number_a, d%dpg_a, d%sigmag_a, nsect_a)
        call size_distribution(1, 1, 1, size(d%number_b), in2b, fn2b, d%number_b, d%dpg_b, d%sigmag_b, nsect_b)
        ndist = nsect_a(1,1,:) + nsect_b(1,1,:)

        a_naerop%d(kp,ip,jp,1:nbins) = ndist

        iso4 = 0
        ioc = 0
        if (spec%isUsed('SO4')) iso4 = spec%getIndex('SO4')
        if (spec%isUsed('OC')) ioc = spec%getIndex('OC')
        if (iso4 > 0 .and. ioc > 0) then
            frac_oc_1a = vf_a(ioc)/(vf_a(ioc) + vf_a(iso4))
        else if (ioc > 0) then
            frac_oc_1a = 1.
        else if (iso4 > 0) then
            frac_oc_1a = 0.
        else
            error stop "salsa_0d: either OC or SO4 must be active for aerosol regime 1a"
        end if

        if (iso4 > 0) a_maerop%d(kp,ip,jp,getMassIndex(nbins,in1a,iso4):getMassIndex(nbins,fn1a,iso4)) = &
            max(0., 1. - frac_oc_1a)*ndist(in1a:fn1a)*core(in1a:fn1a)*spec%rhosu
        if (ioc > 0) a_maerop%d(kp,ip,jp,getMassIndex(nbins,in1a,ioc):getMassIndex(nbins,fn1a,ioc)) = &
            max(0., frac_oc_1a)*ndist(in1a:fn1a)*core(in1a:fn1a)*spec%rhooc

        do s = 1, ndry
            a_maerop%d(kp,ip,jp,getMassIndex(nbins,in2a,s):getMassIndex(nbins,fn2a,s)) = &
                max(0., vf_a(s))*ndist(in2a:fn2a)*core(in2a:fn2a)*spec%rholiq(s)
            a_maerop%d(kp,ip,jp,getMassIndex(nbins,in2b,s):getMassIndex(nbins,fn2b,s)) = &
                max(0., vf_b(s))*ndist(in2b:fn2b)*core(in2b:fn2b)*spec%rholiq(s)
        end do
    end subroutine set_initial_aerosol


    ! Gas tracers, in the SALSA order H2SO4, HNO3, NH3, OCNV, OCSV.
    subroutine set_initial_gases(d)
        type(Input_data), intent(in) :: d
        character(len=8), parameter :: names(5) = [character(len=8) :: 'H2SO4', 'HNO3', 'NH3', 'OCNV', 'OCSV']
        if (size(d%gas) /= 5) error stop "salsa_0d: expected 5 gases in the input file"
        if (any(d%gas_name /= names)) error stop "salsa_0d: unexpected gas order in the input file"
        a_gaerop%d(kp,ip,jp,1:5) = d%gas
    end subroutine set_initial_gases


    function read_text_file(filename) result(text)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: text
        character(len=1024) :: line
        integer :: u, ios
        text = ''
        open(newunit=u, file=filename, status='old', action='read')
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            text = text // trim(line) // achar(10)
        end do
        close(u)
    end function read_text_file

end program salsa_0d
