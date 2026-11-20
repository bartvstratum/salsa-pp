/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Bins, species and gases, and the order of their fields.
 */

#include "salsa_layout.h"


namespace Salsa
{
    Layout::Layout(
            const std::array<int, 2>& n_bins_regime,
            const int n_bins_rain,
            const std::vector<std::string>& dry_species)
    {
        const int n1 = n_bins_regime[0];
        const int n2 = n_bins_regime[1];

        aer_1a_start = 0;
        aer_1a_end = aer_1a_start + n1;
        aer_2a_start = aer_1a_end;
        aer_2a_end = aer_2a_start + n2;
        aer_2b_start = aer_2a_end;
        aer_2b_end = aer_2b_start + n2;
        cld_a_start = aer_2b_end;
        cld_a_end = cld_a_start + n2;
        cld_b_start = cld_a_end;
        cld_b_end = cld_b_start + n2;
        rain_start = cld_b_end;
        rain_end = rain_start + n_bins_rain;
        n_bins = rain_end;

        species = dry_species;
        species.push_back("H2O");
        n_species = static_cast<int>(species.size());

        n_fields = n_bins * (1 + n_species) + n_gases;
    }
}
