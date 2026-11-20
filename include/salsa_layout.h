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

#include "salsa_types.h"


namespace Salsa
{
    // Field order: per bin the number, then the mass of each species; the gases after all bins.
    KOKKOS_INLINE_FUNCTION
    int number_field(const int bin, const int n_species) { return bin*(1 + n_species); }

    KOKKOS_INLINE_FUNCTION
    int mass_field(const int bin, const int species, const int n_species) { return bin*(1 + n_species) + 1 + species; }

    KOKKOS_INLINE_FUNCTION
    int gas_field(const int gas, const int n_bins, const int n_species) { return n_bins*(1 + n_species) + gas; }


    // The bins are ordered as allSALSA in UCLALES-SALSA: aerosol (regimes 1a, 2a, 2b), cloud (a, b), rain.
    // Each bin has a number field, followed by one mass field per species (the dry species, then water).
    // The gas fields come after all bins, in the UCLALES-SALSA order H2SO4, HNO3, NH3, OCNV, OCSV.
    // Bin ranges are [start, end), as kstart and kend in MicroHH.
    class Layout
    {
        public:
            Layout(
                    const std::array<TF, 4>& regime_limits,
                    const std::array<int, 2>& n_bins_regime,
                    const int n_bins_rain,
                    const TF rain_diameter_low,
                    const TF rain_volume_ratio,
                    const std::vector<std::string>& dry_species);

            Layout(const Layout&) = delete;
            Layout& operator=(const Layout&) = delete;
            Layout(Layout&&) = delete;
            Layout& operator=(Layout&&) = delete;

            int number_field(const int bin) const { return Salsa::number_field(bin, n_species); }
            int mass_field(const int bin, const int species) const { return Salsa::mass_field(bin, species, n_species); }
            int gas_field(const int gas) const { return Salsa::gas_field(gas, n_bins, n_species); }

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

            // Per bin: dry volume limits for aerosol and cloud, drop volume limits for rain (m3), and the centre diameter (m).
            std::vector<TF> volume_low;                     // Fortran: vlolim
            std::vector<TF> volume_high;                    // Fortran: vhilim
            std::vector<TF> diameter_mid;                   // Fortran: dmid
            std::vector<TF> number_limit;                   // Fortran: nlim

            // Per species.
            std::vector<std::string> species;
            std::vector<TF> density;                        // Fortran: rholiq

            int n_bins;
            int n_species;
            int n_fields;
    };
}
