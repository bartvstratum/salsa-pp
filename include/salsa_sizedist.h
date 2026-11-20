/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Number per bin of lognormal modes, for initial aerosol.
 */

#pragma once

#include <vector>

#include "salsa_types.h"
#include "salsa_layout.h"


namespace Salsa
{
    // For the aerosol bins [bin_start, bin_end); zero for the other aerosol bins. The unit is that of number.
    std::vector<TF> lognormal_distribution(
            const std::vector<TF>& number,                  // Number of each mode
            const std::vector<TF>& diameter,                // Geometric mean diameter of each mode (m)
            const std::vector<TF>& sigma,                   // Geometric standard deviation of each mode
            const Layout& layout,
            const int bin_start,
            const int bin_end);
}
