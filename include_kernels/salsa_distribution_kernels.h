/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_update.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Update of the size distribution: particles that left their bin move to the neighbouring bin.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"


namespace Salsa
{
    // Fortran: mo_salsa_update.f90: distr_update
    // Moves the particles that have grown or shrunk out of their bin to the neighbouring bin, repeated until the bins
    // are within their limits. Aerosol and cloud bins follow the dry volume, rain bins the wet volume. Within a bin,
    // the particles are spread uniformly over a volume range with the ratio of the bin limits, around the mean volume.
    KOKKOS_INLINE_FUNCTION
    void update_distribution(
            Bin_array& number,                              // m-3
            Bin_species_array& volume,                      // m3 m-3
            const Array_1d<TF>& volume_low,
            const Array_1d<TF>& volume_high,
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& number_limit,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_a_end,
            const int cld_b_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end)
    {
        const int n_dry = n_species - 1;
        const int max_iterations = 100;
        const TF fraction_max = 0.99_tf;
        const TF volume_min = 1.e-30_tf;            // Minimum dry volume of a cloud bin (m3 m-3)

        auto volume_dry = [&](const int b)
        {
            TF v = 0._tf;
            for (int s=0; s<n_dry; ++s)
                v += volume[b][s];
            return v;
        };

        auto volume_wet = [&](const int b)
        {
            TF v = 0._tf;
            for (int s=0; s<n_species; ++s)
                v += volume[b][s];
            return v;
        };

        // Moves a fraction of the volume of all species and a fraction of the number between two bins.
        auto move = [&](const int from, const int to, const TF volume_fraction, const TF number_fraction)
        {
            for (int s=0; s<n_species; ++s)
            {
                volume[to][s] = volume[to][s] + volume_fraction*volume[from][s];
                volume[from][s] = volume[from][s] - volume_fraction*volume[from][s];
            }

            number[to] = number[to] + number_fraction*number[from];
            number[from] = number[from] - number_fraction*number[from];
        };

        // Aerosol, from large to small. The largest 2a bin does not grow into 2b, and the smallest 2b bin shrinks
        // into the largest 1a bin.
        bool within_bins = false;
        int iteration = 0;

        while (!within_bins)
        {
            within_bins = true;

            for (int b=aer_2b_end-2; b>=aer_1a_start; --b)
            {
                if (number[b] > number_limit(b))
                {
                    if (b == aer_2a_end-1)
                        continue;

                    const TF v_dry = volume_dry(b);
                    const TF v_part = v_dry/number[b];

                    if (volume_low(b) > v_part && b == aer_1a_start)
                        continue;

                    // Shrunk: the whole bin moves down.
                    if (volume_low(b) > v_part)
                    {
                        const int b_to = (b == aer_2b_start) ? aer_1a_end-1 : b-1;

                        number[b_to] = number[b_to] + number[b];
                        number[b] = 0._tf;

                        for (int s=0; s<n_species; ++s)
                        {
                            volume[b_to][s] = volume[b_to][s] + volume[b][s];
                            volume[b][s] = 0._tf;
                        }

                        continue;
                    }

                    const TF dm = diameter_mid(b);
                    if (v_part <= pi6*(dm*dm*dm))
                        continue;

                    // Grown: the part above the upper bin limit moves up.
                    const TF v_ratio = volume_high(b)/volume_low(b);
                    const TF v_part_low = 2._tf*v_part/(1._tf + v_ratio);
                    const TF v_part_high = v_ratio*v_part_low;
                    const TF v_excess = 0.5_tf*(v_part_high + volume_high(b));

                    const TF number_fraction = Kokkos::min(1._tf, (v_part_high - volume_high(b))/(v_part_high - v_part_low));

                    // Not grown past the upper limit: with rounding the fraction can be slightly negative just above
                    // pi6*dm^3.
                    if (number_fraction < 0._tf)
                        continue;

                    const TF volume_fraction = Kokkos::min(fraction_max, number_fraction*v_excess/v_part);

                    if (volume_fraction < 0._tf)
                        Kokkos::abort("Error aerosol 0");

                    for (int s=0; s<n_species; ++s)
                    {
                        const TF v_move = number_fraction*number[b]*v_excess*volume[b][s]/v_dry;
                        volume[b+1][s] = volume[b+1][s] + v_move;
                        volume[b][s] = volume[b][s] - v_move;
                    }

                    number[b+1] = number[b+1] + number_fraction*number[b];
                    number[b] = number[b]*(1._tf - number_fraction);
                }

                if (number[b] > number_limit(b))
                {
                    const TF v_part = volume_dry(b)/number[b];
                    within_bins = (volume_low(b) < v_part && v_part < volume_high(b));
                }
            }

            ++iteration;
            if (iteration > max_iterations)
                Kokkos::abort("Aerosol bin update not converged");
        }

        // Cloud, from large to small; the a and b regimes stay apart.
        within_bins = false;
        iteration = 0;

        while (!within_bins)
        {
            within_bins = true;

            for (int b=cld_b_end-1; b>=cld_a_start; --b)
            {
                if (number[b] > number_limit(b) && volume_dry(b) > volume_min)
                {
                    const TF v_part = volume_dry(b)/number[b];
                    const TF v_ratio = volume_high(b)/volume_low(b);
                    const TF v_part_low = 2._tf*v_part/(1._tf + v_ratio);
                    const TF v_part_high = v_ratio*v_part_low;

                    TF v_excess;
                    TF number_fraction;
                    int b_to;

                    if (v_part < volume_low(b) && b != cld_a_start && b != cld_b_start)
                    {
                        v_excess = 0.5_tf*(v_part_low + volume_low(b));
                        number_fraction = Kokkos::min(1._tf, (volume_low(b) - v_part_low)/(v_part_high - v_part_low));
                        b_to = b-1;
                    }
                    else if (v_part > volume_high(b) && b != cld_a_end-1 && b != cld_b_end-1)
                    {
                        v_excess = 0.5_tf*(v_part_high + volume_high(b));
                        number_fraction = Kokkos::min(1._tf, (v_part_high - volume_high(b))/(v_part_high - v_part_low));
                        b_to = b+1;
                    }
                    else
                        continue;

                    const TF volume_fraction = Kokkos::min(fraction_max, number_fraction*v_excess/v_part);

                    if (volume_fraction < 0._tf)
                        Kokkos::abort("Error cloud 0");

                    move(b, b_to, volume_fraction, number_fraction);
                }

                if (number[b] > number_limit(b) && volume_dry(b) > volume_min)
                {
                    const TF v_part = volume_dry(b)/number[b];
                    within_bins = (volume_low(b) < v_part && v_part < volume_high(b));
                }
            }

            ++iteration;
            if (iteration > max_iterations)
                Kokkos::abort("Cloud bin update not converged");
        }

        // Rain, from large to small, following the wet volume.
        within_bins = false;
        iteration = 0;

        while (!within_bins)
        {
            within_bins = true;

            for (int b=rain_end-1; b>=rain_start; --b)
            {
                if (number[b] > number_limit(b))
                {
                    const TF v_part = volume_wet(b)/number[b];
                    const TF v_ratio = volume_high(b)/volume_low(b);
                    const TF v_part_low = 2._tf*v_part/(1._tf + v_ratio);
                    const TF v_part_high = v_ratio*v_part_low;

                    const TF dm = diameter_mid(b);
                    const TF v_mid = pi6*(dm*dm*dm);
                    const TF v_shrunk = v_mid - 0.2_tf*(v_mid - volume_low(b));

                    TF v_excess;
                    TF number_fraction;
                    int b_to;

                    if (v_part < v_shrunk && b != rain_start)
                    {
                        v_excess = 0.5_tf*(v_part_low + volume_low(b));
                        number_fraction = Kokkos::min(1._tf, (volume_low(b) - v_part_low)/(v_part_high - v_part_low));
                        b_to = b-1;
                    }
                    else if (v_part > v_mid && b != rain_end-1)
                    {
                        v_excess = 0.5_tf*(v_part_high + volume_high(b));
                        number_fraction = Kokkos::min(fraction_max, (v_part_high - volume_high(b))/(v_part_high - v_part_low));
                        b_to = b+1;

                        // As for aerosol: with rounding the fraction can be slightly negative just above v_mid.
                        if (number_fraction < 0._tf)
                            continue;
                    }
                    else
                        continue;

                    const TF volume_fraction = Kokkos::min(fraction_max, number_fraction*v_excess/v_part);

                    if (volume_fraction < 0._tf)
                        Kokkos::abort("Error, volc precp 0");

                    move(b, b_to, volume_fraction, number_fraction);
                }

                if (number[b] > number_limit(b))
                {
                    const TF v_part = volume_wet(b)/number[b];
                    within_bins = (volume_low(b) < v_part && v_part < volume_high(b));
                }
            }

            ++iteration;
            if (iteration > max_iterations)
                Kokkos::abort("Precip bin update not converged");
        }
    }
}
