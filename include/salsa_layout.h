/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Bins, species and gases, and the order of their fields.
 */

#pragma once

#include <array>
#include <string>
#include <vector>


namespace Salsa
{
    // The bins are ordered as allSALSA in UCLALES-SALSA: aerosol (regimes 1a, 2a, 2b), cloud (a, b), rain.
    // Each bin has a number field, followed by one mass field per species (the dry species, then water).
    // The gas fields come after all bins, in the UCLALES-SALSA order H2SO4, HNO3, NH3, OCNV, OCSV.
    // Bin ranges are [start, end), as kstart and kend in MicroHH.
    class Layout
    {
        public:
            Layout(
                    const std::array<int, 2>& n_bins_regime,
                    const int n_bins_rain,
                    const std::vector<std::string>& dry_species);

            Layout(const Layout&) = delete;
            Layout& operator=(const Layout&) = delete;
            Layout(Layout&&) = delete;
            Layout& operator=(Layout&&) = delete;

            int number_field(const int bin) const { return bin * (1 + n_species); }
            int mass_field(const int bin, const int species) const { return bin * (1 + n_species) + 1 + species; }
            int gas_field(const int gas) const { return n_bins * (1 + n_species) + gas; }

            static constexpr int n_gases = 5;

            int aer_1a_start;
            int aer_1a_end;
            int aer_2a_start;
            int aer_2a_end;
            int aer_2b_start;
            int aer_2b_end;
            int cld_a_start;
            int cld_a_end;
            int cld_b_start;
            int cld_b_end;
            int rain_start;
            int rain_end;

            std::vector<std::string> species;

            int n_bins;
            int n_species;
            int n_fields;
    };
}
