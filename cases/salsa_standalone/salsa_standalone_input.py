# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Creates the NetCDF input (forcing + initial conditions) of the SALSA-standalone reference cases A to D.
# The forcing comes from input.dat of the original SALSA-standalone (commit 41da7aa) for all cases.

import datetime
import pathlib

import numpy as np
from scipy.io import netcdf_file

case_dir = pathlib.Path(__file__).resolve().parent
input_dat = case_dir.parents[1] / 'fortran' / 'salsa-standalone-reference' / 'input.dat'

# Constants as in SALSA (mo_submctl.f90).
mair = 28.97e-3   # Molar mass of air (kg mol-1).
rg = 8.314        # Molar gas constant (J mol-1 K-1).
avog = 6.0221e23  # Avogadro number (mol-1).

# The old driver reads one row of input.dat per 1 s time step.
dt_row = 1.

# Number of lognormal modes per regime, as in the &salsa namelist.
nmod = 7

dry_species = ['SO4', 'OC', 'BC', 'DU', 'SS', 'NO', 'NH']
gas_names = ['H2SO4', 'HNO3', 'NH3', 'OCNV', 'OCSV']


def char_array(names, strlen):
    return np.array([list(n.ljust(strlen)) for n in names], dtype='S1')


def write_case(case, forcing, modes_a, modes_b, vol_frac_a, vol_frac_b, gas_m3):
    temperature, pressure, rh = forcing
    nt = temperature.size

    # Number concentrations are given per m3 at the start and stored per kg of air,
    # with the dry-air density of the old driver: rho = M_air p / (R T).
    rho0 = mair * pressure[0] / (rg * temperature[0])

    path = case_dir / f'case_{case.lower()}' / 'input.nc'
    path.parent.mkdir(exist_ok=True)

    with netcdf_file(path, 'w', version=2) as nc:
        nc.title = f'SALSA-standalone reference, case {case}: forcing and initial conditions'
        nc.case_name = case
        nc.source = 'forcing: input.dat of https://github.com/UCLALES-SALSA/SALSA-standalone (commit 41da7aa)'
        nc.history = f'{datetime.datetime.now().isoformat(timespec="seconds")}: created by salsa_standalone_input.py'
        nc.Conventions = 'CF-1.8'

        nc.createDimension('time', nt)
        nc.createDimension('mode', nmod)
        nc.createDimension('species_dry', len(dry_species))
        nc.createDimension('gas', len(gas_names))
        nc.createDimension('name_strlen', 8)

        def var(name, dims, data, units, long_name, dtype='d', **attrs):
            v = nc.createVariable(name, dtype, dims)
            v[:] = data
            v.units = units
            v.long_name = long_name
            for key, value in attrs.items():
                setattr(v, key, value)

        var('time', ('time',), np.arange(nt) * dt_row, 's', 'time since start')
        var('temperature', ('time',), temperature, 'K', 'air temperature', standard_name='air_temperature')
        var('pressure', ('time',), pressure, 'Pa', 'air pressure', standard_name='air_pressure')
        var('relative_humidity', ('time',), rh, '1', 'relative humidity (only the first value is used)',
            standard_name='relative_humidity')

        for regime, modes, vol_frac in (('a', modes_a, vol_frac_a), ('b', modes_b, vol_frac_b)):
            n_cm3, dpg_um, sigmag = modes
            label = 'regimes 1a and 2a' if regime == 'a' else 'regime 2b'
            var(f'number_{regime}', ('mode',), n_cm3 * 1e6 / rho0, 'kg-1',
                f'number per kg of air of the initial lognormal modes, {label}')
            var(f'dpg_{regime}', ('mode',), dpg_um * 1e-6, 'm',
                f'geometric mean diameter of the initial lognormal modes, {label}')
            var(f'sigmag_{regime}', ('mode',), sigmag, '1',
                f'geometric standard deviation of the initial lognormal modes, {label}')
            var(f'vol_frac_{regime}', ('species_dry',), vol_frac, '1',
                f'initial dry volume fraction per compound, {label}')
            if regime == 'a':
                nc.variables['vol_frac_a'].description = \
                    'In regime 1a only SO4 and OC are used, in the ratio of their fractions (as UCLALES-SALSA).'

        var('species_dry_name', ('species_dry', 'name_strlen'), char_array(dry_species, 8), '1',
            'dry compound names', dtype='c')
        var('gas_name', ('gas', 'name_strlen'), char_array(gas_names, 8), '1', 'gas names', dtype='c')
        var('gas', ('gas',), gas_m3 / rho0, 'kg-1', 'initial gas concentration, molecules per kg of air')


def main():
    data = np.loadtxt(input_dat)
    forcing = (data[:, 0], data[:, 1] * 100., data[:, 2])

    # driver.f90: two sulphate modes, N = (100, 1) cm-3, Dg = (0.15, 1.0) um, sigma_g = (2.0, 1.5);
    # the other five modes are empty.
    modes_a = (np.array([100., 1., 0., 0., 0., 0., 0.]),
               np.array([0.15, 1.0, 0.2, 0.2, 6.0, 6.0, 6.0]),
               np.array([2.0, 1.5, 1.59, 1.59, 2.0, 2.0, 2.0]))

    # driver.f90: H2SO4 and organic vapours zero, HNO3 = 8.4682e-15 mol cm-3. The old driver also puts the
    # HNO3 value into the NH3 slot; on master neither gas is used.
    c_hno3 = 8.4682e-15 * avog * 1e6
    gas_m3 = np.array([0., c_hno3, c_hno3, 0., 0.])

    # Cases A and B: the initial aerosol of driver.f90, regime 2b empty.
    no_modes = (np.zeros(nmod), modes_a[1].copy(), modes_a[2].copy())
    sulphate = np.array([1., 0., 0., 0., 0., 0., 0.])
    for case in ('A', 'B'):
        write_case(case, forcing, modes_a, no_modes, sulphate, np.zeros(len(dry_species)), gas_m3)

    # Case C: all dry compounds. Mostly soluble a-bins, mostly insoluble (BC, DU) b-bins with one mode of
    # 10 cm-3; its Dg and sigma_g are those of the first b-mode in UCLALES-SALSA bin/runles_dcc.
    modes_b = (np.array([10., 0., 0., 0., 0., 0., 0.]),
               np.array([0.5, 1.0, 0.2, 0.2, 6.0, 6.0, 6.0]),
               np.array([1.7, 1.5, 1.59, 1.59, 2.0, 2.0, 2.0]))
    vol_frac_a = np.array([0.30, 0.25, 0.05, 0.05, 0.15, 0.10, 0.10])
    vol_frac_b = np.array([0.05, 0.05, 0.40, 0.40, 0.04, 0.03, 0.03])
    write_case('C', forcing, modes_a, modes_b, vol_frac_a, vol_frac_b, gas_m3)

    # Case D: case C going up and down, so that rain forms and cloud and rain evaporate back to aerosol.
    # input.dat rows 0-3000 (up to 2000 s, then constant), then rows 2999-0 back down. RH at the start is
    # 95 % instead of 98 %, so that the end is dry enough for the rain to evaporate back to aerosol.
    idx = np.r_[0:3001, 2999:-1:-1]
    rh = forcing[2][idx]
    rh[0] = 0.95
    write_case('D', (forcing[0][idx], forcing[1][idx], rh), modes_a, modes_b, vol_frac_a, vol_frac_b, gas_m3)


if __name__ == '__main__':
    main()
