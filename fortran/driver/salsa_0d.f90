! SALSA++
! Copyright (c) 2026 Wageningen University & Research
! SPDX-License-Identifier: Apache-2.0
!
! 0D driver for the Fortran reference. It uses upstream's own LES glue: the state is set up with
! setPrognosticVariables/setDiagnosticVariables on a 5 x 5 x 3 grid (one interior point),
! run_SALSA computes the SALSA tendencies, and SALSA_diagnostics moves evaporated cloud/rain back to aerosol.
!
! Current status: smoke test. One time step on an empty particle state with fixed thermodynamics;
! forcing, initial aerosol and NetCDF output follow.
program salsa_0d
    use grid, only: nxp, nyp, nzp, level, dtlt
    use classFieldArray, only: FieldArray
    use mo_submctl, only: spec, nbins, ncld, nprc
    use mo_salsa_init, only: define_salsa, salsa_initialize
    use mo_progn_state, only: setPrognosticVariables, a_rp
    use mo_diag_state, only: setDiagnosticVariables, a_temp, a_press, a_dn, a_rsl, a_rh
    use mo_salsa_driver, only: run_SALSA
    use constrain_SALSA, only: SALSA_diagnostics
    implicit none

    real, parameter :: rd = 287.04  ! Gas constant of dry air (J kg-1 K-1).

    type(FieldArray) :: prog, diag
    real, allocatable, target :: a_sclrp(:,:,:,:), a_sclrt(:,:,:,:)
    real, allocatable :: wp(:,:,:), nactd(:,:,:,:), vactd(:,:,:,:)
    character(len=8) :: outputlist(1) = ['none']
    integer :: nwet, nscl, memsize

    ! Reads the &salsa namelist from the file NAMELIST in the working directory.
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
    wp = 0.

    prog = FieldArray()
    diag = FieldArray()
    memsize = 0
    call setPrognosticVariables(a_sclrp, a_sclrt, prog, outputlist, level, 1, .false., nzp, nxp, nyp, nscl)
    call setDiagnosticVariables(diag, outputlist, memsize, level, 0, .false., nzp, nxp, nyp)

    ! Smoke-test thermodynamics: the first row of input.dat.
    a_temp%d = 284.998
    a_press%d = 88637.4
    a_dn%d = a_press%d / (rd*a_temp%d)
    a_rsl%d = 0.622*exp(20.386 - 5132./a_temp%d)*133.32/a_press%d
    a_rh%d = 0.98
    a_rp%d = a_rh%d*a_rsl%d

    call run_SALSA(diag, prog, nzp, nxp, nyp, nwet, wp, nactd, vactd, dtlt, 0., 1, level, .false., .false.)
    a_sclrp = a_sclrp + dtlt*a_sclrt
    a_sclrt = 0.
    call SALSA_diagnostics(.false., .false.)

    print '(a,i0,a,i0,a,i0,a,i0,a)', 'salsa_0d smoke test OK: ', nbins, ' aerosol, ', ncld, ' cloud, ', &
        nprc, ' rain bins, ', nwet, ' wet compounds'
end program salsa_0d
