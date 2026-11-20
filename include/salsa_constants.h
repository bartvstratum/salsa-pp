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

    constexpr TF nlim = 1._tf;              // Number limit for aerosol and cloud (kg-1)
    constexpr TF prlim = 1.e-6_tf;          // Number limit for rain (kg-1)
}
