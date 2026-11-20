! SALSA++
! Copyright (c) 2026 Wageningen University & Research
! SPDX-License-Identifier: Apache-2.0
!
! NetCDF input and output of the 0D driver. The output layout is shared with the C++ driver.
module salsa_0d_io
    use netcdf
    implicit none
    private
    public :: Input_data, read_input, create_output, write_output, close_output

    ! Grid point of the 5 x 5 x 3 grid where SALSA runs.
    integer, parameter :: kp = 2, ip = 3, jp = 3

    integer, parameter :: name_strlen = 8

    type Input_data
        real, allocatable :: time(:), temperature(:), pressure(:), relative_humidity(:)
        real, allocatable :: number_a(:), dpg_a(:), sigmag_a(:), number_b(:), dpg_b(:), sigmag_b(:)
        character(len=name_strlen), allocatable :: species_dry_name(:), gas_name(:)
        real, allocatable :: vol_frac_a(:), vol_frac_b(:), gas(:)
    end type Input_data

    ! Output file state.
    integer :: ncid = -1
    integer :: irec = 0
    integer :: varid_time, varid_temperature, varid_pressure, varid_density, varid_qv, varid_qsat, varid_ss, varid_gas
    integer :: varid_number(3), varid_mass(3), varid_dwet(3), varid_ddry(3)
    character(len=*), parameter :: category(3) = [character(len=7) :: 'aerosol', 'cloud', 'rain']

contains

    subroutine check(status, what)
        integer, intent(in) :: status
        character(len=*), intent(in) :: what
        if (status /= nf90_noerr) then
            print '(a)', 'NetCDF error (' // what // '): ' // trim(nf90_strerror(status))
            error stop 1
        end if
    end subroutine check


    subroutine read_input(filename, d)
        character(len=*), intent(in) :: filename
        type(Input_data), intent(out) :: d
        integer :: id

        call check(nf90_open(filename, nf90_nowrite, id), filename)
        call read_1d(id, 'time', d%time)
        call read_1d(id, 'temperature', d%temperature)
        call read_1d(id, 'pressure', d%pressure)
        call read_1d(id, 'relative_humidity', d%relative_humidity)
        call read_1d(id, 'number_a', d%number_a)
        call read_1d(id, 'dpg_a', d%dpg_a)
        call read_1d(id, 'sigmag_a', d%sigmag_a)
        call read_1d(id, 'number_b', d%number_b)
        call read_1d(id, 'dpg_b', d%dpg_b)
        call read_1d(id, 'sigmag_b', d%sigmag_b)
        call read_1d(id, 'vol_frac_a', d%vol_frac_a)
        call read_1d(id, 'vol_frac_b', d%vol_frac_b)
        call read_1d(id, 'gas', d%gas)
        call read_names(id, 'species_dry_name', d%species_dry_name)
        call read_names(id, 'gas_name', d%gas_name)
        call check(nf90_close(id), filename)
    end subroutine read_input


    subroutine read_1d(id, name, values)
        integer, intent(in) :: id
        character(len=*), intent(in) :: name
        real, allocatable, intent(out) :: values(:)
        integer :: varid, dimids(1), n

        call check(nf90_inq_varid(id, name, varid), name)
        call check(nf90_inquire_variable(id, varid, dimids=dimids), name)
        call check(nf90_inquire_dimension(id, dimids(1), len=n), name)
        allocate(values(n))
        call check(nf90_get_var(id, varid, values), name)
    end subroutine read_1d


    subroutine read_names(id, name, names)
        integer, intent(in) :: id
        character(len=*), intent(in) :: name
        character(len=name_strlen), allocatable, intent(out) :: names(:)
        integer :: varid, dimids(2), n

        call check(nf90_inq_varid(id, name, varid), name)
        call check(nf90_inquire_variable(id, varid, dimids=dimids), name)
        call check(nf90_inquire_dimension(id, dimids(2), len=n), name)
        allocate(names(n))
        call check(nf90_get_var(id, varid, names), name)
    end subroutine read_names


    ! Creates the output file and writes the coordinates. Settings go into global attributes.
    subroutine create_output(filename, case_name, input_file, namelist_text, git_commit)
        use mo_submctl, only: spec, nbins, ncld, nprc, pi6
        use mo_salsa_types, only: aero, cloud, precp
        character(len=*), intent(in) :: filename, case_name, input_file, namelist_text, git_commit

        integer :: dimid_time, dimid_bin(3), dimid_species, dimid_gas, dimid_nv, dimid_strlen
        integer :: varid_bin(3), varid_bounds(3), varid_species, varid_gas_name
        integer :: nbin(3), nwet, ic, nb
        character(len=19) :: creation_time
        integer :: dt_values(8)
        real, allocatable :: dmid(:), bounds(:,:)
        character(len=name_strlen), allocatable :: names(:)
        character(len=*), parameter :: gas_names(5) = &
            [character(len=name_strlen) :: 'H2SO4', 'HNO3', 'NH3', 'OCNV', 'OCSV']

        nbin = [nbins, ncld, nprc]
        nwet = spec%getNSpec(type="wet")

        call check(nf90_create(filename, nf90_64bit_offset, ncid), filename)

        call date_and_time(values=dt_values)
        write(creation_time, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') dt_values(1:3), dt_values(5:7)

        call put_att_global('title', 'SALSA 0D reference test, case ' // trim(case_name))
        call put_att_global('case_name', trim(case_name))
        call put_att_global('code', 'Fortran reference: UCLALES-SALSA (vendored) + SALSA++ Fortran driver')
        call put_att_global('salsa_pp_git_commit', trim(git_commit))
        call put_att_global('upstream_salsa_commit', &
            'https://github.com/UCLALES-SALSA/UCLALES-SALSA 8c779581118f97b6092efe6da43680ffa4e1cefa')
        call put_att_global('creation_time', creation_time)
        call put_att_global('input_file', trim(input_file))
        call put_att_global('namelist', namelist_text)
        call put_att_global('Conventions', 'CF-1.8')

        call check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_time), 'time')
        do ic = 1, 3
            call check(nf90_def_dim(ncid, 'bin_' // trim(category(ic)), nbin(ic), dimid_bin(ic)), 'bin')
        end do
        call check(nf90_def_dim(ncid, 'species', nwet, dimid_species), 'species')
        call check(nf90_def_dim(ncid, 'gas', 5, dimid_gas), 'gas')
        call check(nf90_def_dim(ncid, 'nv', 2, dimid_nv), 'nv')
        call check(nf90_def_dim(ncid, 'name_strlen', name_strlen, dimid_strlen), 'name_strlen')

        varid_time = def_var('time', [dimid_time], 's', 'time since start')

        ! Bin coordinates: dry diameter for aerosol and cloud (cloud bins are parallel to aerosol bins), drop diameter for rain.
        do ic = 1, 3
            if (ic < 3) then
                varid_bin(ic) = def_var('bin_' // trim(category(ic)), [dimid_bin(ic)], 'm', &
                    trim(category(ic)) // ' bin centre dry diameter (SALSA dmid)')
            else
                varid_bin(ic) = def_var('bin_rain', [dimid_bin(ic)], 'm', 'rain bin centre drop diameter (SALSA dmid)')
            end if
            call check(nf90_put_att(ncid, varid_bin(ic), 'bounds', 'bin_' // trim(category(ic)) // '_bounds'), 'bounds')
            varid_bounds(ic) = def_var('bin_' // trim(category(ic)) // '_bounds', [dimid_nv, dimid_bin(ic)], 'm', &
                trim(category(ic)) // ' bin diameter limits')
        end do

        call check(nf90_def_var(ncid, 'species', nf90_char, [dimid_strlen, dimid_species], varid_species), 'species')
        call check(nf90_put_att(ncid, varid_species, 'long_name', 'compound names, in the order of the species dimension'), &
            'species')
        call check(nf90_def_var(ncid, 'gas_name', nf90_char, [dimid_strlen, dimid_gas], varid_gas_name), 'gas_name')
        call check(nf90_put_att(ncid, varid_gas_name, 'long_name', 'gas names, in the order of the gas dimension'), 'gas_name')

        varid_temperature = def_var('temperature', [dimid_time], 'K', 'air temperature', 'air_temperature')
        varid_pressure = def_var('pressure', [dimid_time], 'Pa', 'air pressure', 'air_pressure')
        varid_density = def_var('air_density', [dimid_time], 'kg m-3', 'dry air density, M_air p / (R T)')
        varid_qv = def_var('qv', [dimid_time], 'kg kg-1', 'water vapour mixing ratio', 'humidity_mixing_ratio')
        varid_qsat = def_var('qsat', [dimid_time], 'kg kg-1', &
            'saturation water vapour mixing ratio, M_w e_s(T) / (M_air p), e_s from a 6th-order polynomial')
        varid_ss = def_var('supersaturation', [dimid_time], '1', 'supersaturation with respect to liquid water, qv/qsat - 1')
        varid_gas = def_var('gas', [dimid_gas, dimid_time], 'kg-1', 'gas concentration, molecules per kg of air')

        do ic = 1, 3
            varid_number(ic) = def_var('number_' // trim(category(ic)), [dimid_bin(ic), dimid_time], 'kg-1', &
                trim(category(ic)) // ' number per kg of air')
            varid_mass(ic) = def_var('mass_' // trim(category(ic)), [dimid_bin(ic), dimid_species, dimid_time], 'kg kg-1', &
                trim(category(ic)) // ' mass per compound per kg of air')
            varid_dwet(ic) = def_var('diameter_wet_' // trim(category(ic)), [dimid_bin(ic), dimid_time], 'm', &
                trim(category(ic)) // ' wet diameter: volume-equivalent diameter of the mean particle, all compounds')
            varid_ddry(ic) = def_var('diameter_dry_' // trim(category(ic)), [dimid_bin(ic), dimid_time], 'm', &
                trim(category(ic)) // ' dry diameter: volume-equivalent diameter of the mean particle, dry compounds')
            call check(nf90_put_att(ncid, varid_dwet(ic), '_FillValue', nf90_fill_double), 'fill')
            call check(nf90_put_att(ncid, varid_ddry(ic), '_FillValue', nf90_fill_double), 'fill')
            call check(nf90_put_att(ncid, varid_dwet(ic), 'description', &
                'Fill value where the bin is empty (number below the SALSA limit nlim, or prlim for rain).'), 'description')
        end do

        call check(nf90_enddef(ncid), 'enddef')

        ! Coordinates.
        do ic = 1, 3
            allocate(dmid(nbin(ic)), bounds(2, nbin(ic)))
            do nb = 1, nbin(ic)
                select case (ic)
                case (1)
                    dmid(nb) = aero(1,1,nb)%dmid
                    bounds(:,nb) = ([aero(1,1,nb)%vlolim, aero(1,1,nb)%vhilim]/pi6)**(1./3.)
                case (2)
                    dmid(nb) = cloud(1,1,nb)%dmid
                    bounds(:,nb) = ([cloud(1,1,nb)%vlolim, cloud(1,1,nb)%vhilim]/pi6)**(1./3.)
                case (3)
                    dmid(nb) = precp(1,1,nb)%dmid
                    bounds(:,nb) = ([precp(1,1,nb)%vlolim, precp(1,1,nb)%vhilim]/pi6)**(1./3.)
                end select
            end do
            call check(nf90_put_var(ncid, varid_bin(ic), dmid), 'bin')
            call check(nf90_put_var(ncid, varid_bounds(ic), bounds), 'bounds')
            deallocate(dmid, bounds)
        end do

        allocate(names(nwet))
        names = spec%names(1:nwet)
        call check(nf90_put_var(ncid, varid_species, names), 'species')
        call check(nf90_put_var(ncid, varid_gas_name, gas_names), 'gas_name')
    end subroutine create_output


    subroutine put_att_global(name, value)
        character(len=*), intent(in) :: name, value
        call check(nf90_put_att(ncid, nf90_global, name, value), name)
    end subroutine put_att_global


    integer function def_var(name, dimids, units, long_name, standard_name)
        character(len=*), intent(in) :: name, units, long_name
        integer, intent(in) :: dimids(:)
        character(len=*), intent(in), optional :: standard_name

        call check(nf90_def_var(ncid, name, nf90_double, dimids, def_var), name)
        call check(nf90_put_att(ncid, def_var, 'units', units), name)
        call check(nf90_put_att(ncid, def_var, 'long_name', long_name), name)
        if (present(standard_name)) &
            call check(nf90_put_att(ncid, def_var, 'standard_name', standard_name), name)
    end function def_var


    ! Writes the state of the SALSA grid point as the next time record.
    subroutine write_output(time)
        use mo_submctl, only: spec, nbins, ncld, nprc, nlim, prlim, pi6
        use mo_progn_state, only: a_rp, a_gaerop, a_naerop, a_maerop, a_ncloudp, a_mcloudp, a_nprecpp, a_mprecpp
        use mo_diag_state, only: a_temp, a_press, a_dn, a_rsl
        real, intent(in) :: time

        integer :: nwet, ndry

        nwet = spec%getNSpec(type="wet")
        ndry = spec%getNSpec(type="dry")
        irec = irec + 1

        call put_scalar(varid_time, time)
        call put_scalar(varid_temperature, a_temp%d(kp,ip,jp))
        call put_scalar(varid_pressure, a_press%d(kp,ip,jp))
        call put_scalar(varid_density, a_dn%d(kp,ip,jp))
        call put_scalar(varid_qv, a_rp%d(kp,ip,jp))
        call put_scalar(varid_qsat, a_rsl%d(kp,ip,jp))
        call put_scalar(varid_ss, a_rp%d(kp,ip,jp)/a_rsl%d(kp,ip,jp) - 1.)
        call check(nf90_put_var(ncid, varid_gas, a_gaerop%d(kp,ip,jp,1:5), start=[1,irec], count=[5,1]), 'gas')

        call put_category(1, nbins, nlim,  a_naerop%d(kp,ip,jp,:), a_maerop%d(kp,ip,jp,:))
        call put_category(2, ncld,  nlim,  a_ncloudp%d(kp,ip,jp,:), a_mcloudp%d(kp,ip,jp,:))
        call put_category(3, nprc,  prlim, a_nprecpp%d(kp,ip,jp,:), a_mprecpp%d(kp,ip,jp,:))

    contains

        subroutine put_scalar(varid, value)
            integer, intent(in) :: varid
            real, intent(in) :: value
            call check(nf90_put_var(ncid, varid, [value], start=[irec], count=[1]), 'scalar')
        end subroutine put_scalar

        ! Masses are stored per compound, with all bins of a compound consecutive.
        subroutine put_category(ic, nb, numlim, number, mass_flat)
            integer, intent(in) :: ic, nb
            real, intent(in) :: numlim, number(:), mass_flat(:)
            real :: mass(nb,nwet), dwet(nb), ddry(nb)
            integer :: b

            mass = reshape(mass_flat(1:nb*nwet), [nb, nwet])
            do b = 1, nb
                if (number(b) > numlim) then
                    dwet(b) = (sum(mass(b,1:nwet)/spec%rholiq(1:nwet))/number(b)/pi6)**(1./3.)
                    ddry(b) = (sum(mass(b,1:ndry)/spec%rholiq(1:ndry))/number(b)/pi6)**(1./3.)
                else
                    dwet(b) = nf90_fill_double
                    ddry(b) = nf90_fill_double
                end if
            end do

            call check(nf90_put_var(ncid, varid_number(ic), number(1:nb), start=[1,irec], count=[nb,1]), 'number')
            call check(nf90_put_var(ncid, varid_mass(ic), mass, start=[1,1,irec], count=[nb,nwet,1]), 'mass')
            call check(nf90_put_var(ncid, varid_dwet(ic), dwet, start=[1,irec], count=[nb,1]), 'dwet')
            call check(nf90_put_var(ncid, varid_ddry(ic), ddry, start=[1,irec], count=[nb,1]), 'ddry')
        end subroutine put_category

    end subroutine write_output


    subroutine close_output()
        call check(nf90_close(ncid), 'close')
        ncid = -1
    end subroutine close_output

end module salsa_0d_io
