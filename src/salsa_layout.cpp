/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_init.f90, classSpecies.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Bins, species and gases, and the order of their fields.
 */

#include <cmath>
#include <map>
#include <stdexcept>

#include "salsa_layout.h"
#include "salsa_constants.h"


namespace
{
    using namespace Salsa;

    // Fortran: classSpecies.f90: allRho
    const std::map<std::string, TF> densities = {
            {"SO4", 1830._tf},
            {"OC", 1400._tf},
            {"BC", 2000._tf},
            {"DU", 2650._tf},
            {"SS", 2165._tf},
            {"NO", 1479._tf},
            {"NH", 1530._tf},
            {"H2O", 1000._tf}};


    // Fortran: classSpecies.f90: allMM
    const std::map<std::string, TF> molar_masses = {
            {"SO4", 98.08e-3_tf},
            {"OC", 150.e-3_tf},
            {"BC", 12.e-3_tf},
            {"DU", 100.e-3_tf},
            {"SS", 58.44e-3_tf},
            {"NO", 62.01e-3_tf},
            {"NH", 18.04e-3_tf},
            {"H2O", 18.016e-3_tf}};


    // Fortran: classSpecies.f90: allDiss
    const std::map<std::string, TF> dissociations = {
            {"SO4", 3._tf},
            {"OC", 1._tf},
            {"BC", 0._tf},
            {"DU", 0._tf},
            {"SS", 2._tf},
            {"NO", 1._tf},
            {"NH", 1._tf},
            {"H2O", 1._tf}};


    // Fortran: mo_salsa_init.f90: set_aerobins (part)
    // Logarithmically spaced bins between two diameters.
    void set_log_bins(
            std::vector<TF>& volume_low,
            std::vector<TF>& volume_high,
            const TF diameter_low,
            const TF diameter_high,
            const int start,
            const int n)
    {
        const TF ratio = diameter_high/diameter_low;

        for (int b=0; b<n; ++b)
        {
            const TF d_low = diameter_low*std::pow(ratio, static_cast<TF>(b)/n);
            const TF d_high = diameter_low*std::pow(ratio, static_cast<TF>(b+1)/n);

            volume_low[start + b] = pi6*(d_low*d_low*d_low);
            volume_high[start + b] = pi6*(d_high*d_high*d_high);
        }
    }


    // Fortran: mo_salsa_init.f90: buildBinLimits
    // Fortran: mo_salsa_init.f90: set_cloudbins (part)
    // Bins with a constant volume ratio, starting at a given diameter.
    void set_geometric_bins(
            std::vector<TF>& volume_low,
            std::vector<TF>& volume_high,
            const TF diameter_low,
            const TF volume_ratio,
            const int start,
            const int n)
    {
        const TF diameter_ratio = std::pow(volume_ratio, 1._tf/3._tf);

        TF d_low = diameter_low;
        for (int b=0; b<n; ++b)
        {
            if (b > 0)
                d_low *= diameter_ratio;

            const TF d_high = d_low*diameter_ratio;

            volume_low[start + b] = pi6*(d_low*d_low*d_low);
            volume_high[start + b] = pi6*(d_high*d_high*d_high);
        }
    }
}


namespace Salsa
{
    Layout::Layout(
            const std::array<TF, 4>& regime_limits,
            const std::array<int, 2>& n_bins_regime,
            const int n_bins_rain,
            const TF rain_diameter_low,
            const TF rain_volume_ratio,
            const std::vector<std::string>& dry_species)
    {
        const int n1 = n_bins_regime[0];
        const int n2 = n_bins_regime[1];

        // UCLALES-SALSA puts the first 4 bins of regime 2 below regime_limits[2].
        const int n2_low = 4;               // Fortran: nbin2
        const int n2_high = n2 - n2_low;    // Fortran: nbin3

        if (n2_high < 0)
            throw std::runtime_error("nbin(2) must be at least 4");

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

        volume_low.resize(n_bins);
        volume_high.resize(n_bins);
        diameter_mid.resize(n_bins);

        set_log_bins(volume_low, volume_high, regime_limits[0], regime_limits[1], aer_1a_start, n1);
        set_log_bins(volume_low, volume_high, regime_limits[1], regime_limits[2], aer_2a_start, n2_low);
        set_log_bins(volume_low, volume_high, regime_limits[2], regime_limits[3], aer_2a_start + n2_low, n2_high);

        // Regime 2b and the cloud bins (a and b) are parallel to regime 2a.
        for (int b=0; b<n2; ++b)
        {
            for (const int start : {aer_2b_start, cld_a_start, cld_b_start})
            {
                volume_low[start + b] = volume_low[aer_2a_start + b];
                volume_high[start + b] = volume_high[aer_2a_start + b];
            }
        }

        set_geometric_bins(volume_low, volume_high, rain_diameter_low, rain_volume_ratio, rain_start, n_bins_rain);

        // Fortran: mo_salsa_init.f90: set_aerobins (part)
        for (int b=0; b<n_bins; ++b)
            diameter_mid[b] = std::pow((volume_high[b] + volume_low[b])/(2._tf*pi6), 1._tf/3._tf);

        // Fortran: mo_salsa_init.f90: set_aerobins (part), set_cloudbins (part)
        number_limit.assign(n_bins, nlim);
        diameter_limit.assign(n_bins, dlaero);
        for (int b=cld_a_start; b<cld_b_end; ++b)
            diameter_limit[b] = dlcloud;
        for (int b=rain_start; b<rain_end; ++b)
        {
            number_limit[b] = prlim;
            diameter_limit[b] = dlprecp;
        }

        species = dry_species;
        species.push_back("H2O");
        n_species = static_cast<int>(species.size());

        for (const std::string& name : species)
        {
            const auto rho = densities.find(name);
            if (rho == densities.end())
                throw std::runtime_error("Unknown species " + name);

            density.push_back(rho->second);
            molar_mass.push_back(molar_masses.at(name));
            dissociation.push_back(dissociations.at(name));
        }

        n_fields = n_bins * (1 + n_species) + Gas::count;

        if (n_bins > max_bins)
            throw std::runtime_error("More bins than max_bins in salsa_types.h");
        if (n_species > max_species)
            throw std::runtime_error("More species than max_species in salsa_types.h");
    }
}
