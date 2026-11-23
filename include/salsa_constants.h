/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_submctl.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Physical constants and number limits, with the values of UCLALES-SALSA (pi and pi6 truncated as there).
 */

#pragma once

#include "salsa_types.h"


namespace Salsa
{
    constexpr TF pi = 3.1415927_tf;
    constexpr TF pi6 = 0.5235988_tf;
    constexpr TF avog = 6.0221e+23_tf;      // Avogadro number (mol-1)
    constexpr TF boltz = 1.3807e-23_tf;     // Boltzmann constant (J K-1)
    constexpr TF rg = 8.314_tf;             // Molar gas constant (J mol-1 K-1)
    constexpr TF mair = 28.97e-3_tf;        // Molar mass of air (kg mol-1)
    constexpr TF mwa = 18.016e-3_tf;        // Molar mass of water (kg mol-1)
    constexpr TF rhowa = 1000._tf;          // Density of water (kg m-3)
    constexpr TF mvwa = mwa/avog/rhowa;     // Volume of one water molecule (m3)
    constexpr TF surfw0 = 0.073_tf;         // Surface tension of water (J m-2)
    constexpr TF epsoc = 0.15_tf;           // Water uptake of organic carbon (fraction)
    constexpr TF alv = 2.5e6_tf;            // Latent heat of vaporisation (J kg-1)
    constexpr TF pstand = 1.01325e+5_tf;    // Standard pressure (Pa)
    constexpr TF grav = 9.81_tf;            // Gravitational acceleration (m s-2)
    constexpr TF rd = 287.04_tf;            // Gas constant of dry air (J kg-1 K-1)
    constexpr TF eddy_dis = 0.001_tf;       // Dissipation rate of turbulent kinetic energy (m2 s-3)
    constexpr TF eps = static_cast<TF>(Kokkos::Experimental::epsilon_v<double>);  // Also in float: used as an absolute threshold

    constexpr TF msu = 98.08e-3_tf;         // Molar mass of sulphate (kg mol-1)
    constexpr TF moc = 150.e-3_tf;          // Molar mass of organic carbon (kg mol-1)
    constexpr TF mno = 62.01e-3_tf;         // Molar mass of HNO3 (kg mol-1)
    constexpr TF mnh = 18.04e-3_tf;         // Molar mass of NH3 (kg mol-1)
    constexpr TF rhosu = 1830._tf;          // Density of sulphate (kg m-3)
    constexpr TF rhooc = 1400._tf;          // Density of organic carbon (kg m-3)
    constexpr TF mvsu = msu/avog/rhosu;     // Volume of one sulphate molecule (m3)
    constexpr TF mvoc = moc/avog/rhooc;     // Volume of one organic carbon molecule (m3)
    constexpr TF d_sa = 5.539376964394570e-10_tf;   // Diameter of a H2SO4 molecule (m)

    constexpr TF dlaero = 30.e-6_tf;        // Maximum wet diameter of aerosol (m)
    constexpr TF dlcloud = 100.e-6_tf;      // Maximum wet diameter of cloud droplets (m)
    constexpr TF dlprecp = 5.e-3_tf;        // Maximum wet diameter of rain drops (m)

    constexpr TF nlim = 1._tf;              // Number limit for aerosol and cloud (kg-1)
    constexpr TF prlim = 1.e-6_tf;          // Number limit for rain (kg-1)
}
