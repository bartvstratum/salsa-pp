/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_cloud.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Activation of aerosol to cloud droplets.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"


namespace Salsa
{
    // Fortran: mo_salsa_cloud.f90: intgN
    // Integral of a linear number distribution (slope, origin) between two dry volumes.
    KOKKOS_INLINE_FUNCTION
    TF integral_number(
            const TF slope,
            const TF origin,
            const TF low,
            const TF high)
    {
        return 0.5_tf*slope*Kokkos::max(high*high - low*low, 0._tf) + origin*Kokkos::max(high - low, 0._tf);
    }


    // Fortran: mo_salsa_cloud.f90: intgV
    // Integral of the volume of a linear number distribution (slope, origin) between two dry volumes.
    KOKKOS_INLINE_FUNCTION
    TF integral_volume(
            const TF slope,
            const TF origin,
            const TF low,
            const TF high)
    {
        return (1._tf/3._tf)*slope*Kokkos::max(high*high*high - low*low*low, 0._tf)
             + 0.5_tf*origin*Kokkos::max(high*high - low*low, 0._tf);
    }


    // Fortran: mo_salsa_cloud.f90: actInterst
    // Activation of the regime 2a and 2b aerosol whose wet size exceeds the critical size at the current
    // supersaturation, into the parallel cloud bins. Within each bin, the wet volume and the number density are
    // linear profiles through the neighbouring bins.
    KOKKOS_INLINE_FUNCTION
    void activate(
            Bin_array& number,                              // m-3
            Bin_species_array& volume,                      // m3 m-3
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& number_limit,
            const Array_1d<TF>& volume_ratio_low,
            const Array_1d<TF>& volume_ratio_high,
            const Array_1d<bool>& soluble,
            const TF T,
            const TF rv,
            const TF rs,
            const int iwa,
            const int n_species,
            const int aer_1a_start,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_a_end,
            const int cld_b_start,
            const int cld_b_end)
    {
        if (rv/rs <= 1._tf)
            return;

        const int n_dry = n_species - 1;
        const TF volume_threshold = 1.e-28_tf;      // Less than the volume of a 1 nm particle (m3)
        const TF none = 999._tf;

        const TF kelvin = 4._tf*mwa*surfw0/(rg*rhowa*T);
        const TF d_crit = 2._tf*kelvin/(3._tf*((rv/rs) - 1._tf));
        const TF v_crit = pi6*(d_crit*d_crit*d_crit);

        auto parallel_aer = [&](const int cb)
        {
            return (cb < cld_a_end) ? aer_2a_start + (cb - cld_a_start) : aer_2b_start + (cb - cld_b_start);
        };

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

        // Activated number and volume per cloud bin; all bins first, then the transfer.
        Bin_array number_act;
        Bin_species_array volume_act;

        for (int cb=cld_a_start; cb<cld_b_end; ++cb)
        {
            const int ab = parallel_aer(cb);

            number_act[cb] = 0._tf;
            for (int s=0; s<n_species; ++s)
                volume_act[cb][s] = 0._tf;

            TF volume_soluble = 0._tf;
            for (int s=0; s<n_species; ++s)
                if (soluble(s))
                    volume_soluble += volume[ab][s];

            if (number[ab] < number_limit(ab))
                continue;
            if (volume_soluble < number[ab]*volume_threshold)
                continue;
            if (volume[ab][iwa] < number[ab]*volume_threshold)
                continue;

            bool range[4] = {false, false, false, false};

            // Dry and wet volume of a particle, in the current bin and at its limits.
            const TF n_mid = number[ab];
            const TF vw_mid = volume_wet(ab)/n_mid;
            const TF v_mid = volume_dry(ab)/n_mid;
            const TF v_lo = v_mid*volume_ratio_low(ab);
            const TF v_hi = v_mid*volume_ratio_high(ab);

            // The previous (smaller) bin, or values for it at the start of a regime.
            TF n_im1;
            TF v_im1;
            TF vw_im1;
            TF v_lo_im1;
            TF v_hi_im1;

            if (ab == aer_1a_start || ab == aer_2b_start)
            {
                n_im1 = number_limit(ab);
                v_im1 = v_lo/2._tf;
                v_lo_im1 = 0._tf;
                v_hi_im1 = v_lo;
                vw_im1 = v_im1;
            }
            else
            {
                n_im1 = number[ab-1];
                if (n_im1 > number_limit(ab-1))
                {
                    v_im1 = volume_dry(ab-1)/n_im1;
                    vw_im1 = volume_wet(ab-1)/n_im1;
                }
                else
                {
                    const TF d = diameter_mid(ab-1);
                    v_im1 = pi6*(d*d*d);
                    vw_im1 = pi6*(d*d*d);
                }
                v_lo_im1 = v_im1*volume_ratio_low(ab-1);
                v_hi_im1 = v_im1*volume_ratio_high(ab-1);
            }

            // The next (larger) bin, or values for it at the end of a regime.
            TF n_ip1;
            TF v_ip1;
            TF vw_ip1;
            TF v_lo_ip1;
            TF v_hi_ip1;

            if (ab == aer_2a_end - 1 || ab == aer_2b_end - 1)
            {
                n_ip1 = number_limit(ab);
                v_ip1 = v_hi + 0.5_tf*(v_hi - v_lo);
                v_lo_ip1 = v_hi;
                v_hi_ip1 = v_hi + (v_hi - v_lo);
                vw_ip1 = v_ip1;
            }
            else
            {
                n_ip1 = number[ab+1];
                if (n_ip1 > number_limit(ab+1))
                {
                    v_ip1 = volume_dry(ab+1)/n_ip1;
                    vw_ip1 = volume_wet(ab+1)/n_ip1;
                }
                else
                {
                    const TF d = diameter_mid(ab+1);
                    v_ip1 = pi6*(d*d*d);
                    vw_ip1 = pi6*(d*d*d);
                }
                v_lo_ip1 = v_ip1*volume_ratio_low(ab+1);
                v_hi_ip1 = v_ip1*volume_ratio_high(ab+1);
            }

            v_ip1 = Kokkos::max(v_hi, v_ip1);
            v_im1 = Kokkos::min(v_lo, v_im1);

            // Wet volume profile: slopes to the neighbouring bin centres, interpolated to the bin limits.
            TF slope_1 = (vw_mid - Kokkos::max(vw_im1, 0._tf))/(v_mid - v_im1);
            TF slope_2 = (Kokkos::max(vw_ip1, 0._tf) - vw_mid)/(v_ip1 - v_mid);

            const TF vw_origin_1 = vw_mid - slope_1*v_mid;
            const TF vw_origin_2 = vw_mid - slope_2*v_mid;

            const TF vw_lo = Kokkos::max(vw_origin_1 + slope_1*v_lo, v_lo);
            const TF vw_hi = Kokkos::max(vw_origin_2 + slope_2*v_hi, v_hi);

            // Dry volume limits of the activating part.
            TF v_int_1;
            TF v_int_2;

            if (v_crit < vw_lo && v_crit < vw_mid && v_crit < vw_hi)
            {
                v_int_1 = v_lo;
                v_int_2 = v_hi;
                for (int r=0; r<4; ++r)
                    range[r] = true;
            }
            else if (v_crit > vw_lo && v_crit > vw_mid && v_crit > vw_hi)
            {
                v_int_1 = none;
                v_int_2 = none;
            }
            else
            {
                TF v_cut = (v_crit - vw_origin_1)/slope_1;
                v_int_1 = (v_cut < v_lo || v_cut > v_mid) ? v_lo : v_cut;

                v_cut = (v_crit - vw_origin_2)/slope_2;
                v_int_2 = (v_cut < v_mid || v_cut > v_hi) ? v_mid : v_cut;

                range[0] = (vw_lo > v_crit);
                range[1] = (vw_mid > v_crit);
                range[2] = (vw_mid > v_crit);
                range[3] = (vw_hi > v_crit);
            }

            // Number density profile and its integrals.
            const TF dn_im1 = n_im1/(v_hi_im1 - v_lo_im1);
            const TF dn_ip1 = n_ip1/(v_hi_ip1 - v_lo_ip1);
            const TF dn_mid = n_mid/(v_hi - v_lo);

            slope_1 = (dn_mid - dn_im1)/(v_mid - v_im1);
            slope_2 = (dn_ip1 - dn_mid)/(v_ip1 - v_mid);

            const TF n_origin_1 = dn_mid - slope_1*v_mid;
            const TF n_origin_2 = dn_mid - slope_2*v_mid;

            const TF n_norm = integral_number(slope_1, n_origin_1, v_lo, v_mid)
                            + integral_number(slope_2, n_origin_2, v_mid, v_hi);
            const TF v_norm = integral_volume(slope_1, n_origin_1, v_lo, v_mid)
                            + integral_volume(slope_2, n_origin_2, v_mid, v_hi);

            TF n_act = 0._tf;
            TF v_act = 0._tf;

            if (range[0])
            {
                n_act = n_act + (n_mid/n_norm)*integral_number(slope_1, n_origin_1, v_lo, v_int_1);
                v_act = v_act + (n_mid*v_mid/v_norm)*integral_volume(slope_1, n_origin_1, v_lo, v_int_1);
            }

            if (range[1])
            {
                n_act = n_act + (n_mid/n_norm)*integral_number(slope_1, n_origin_1, v_int_1, v_mid);
                v_act = v_act + (n_mid*v_mid/v_norm)*integral_volume(slope_1, n_origin_1, v_int_1, v_mid);
            }

            if (range[2])
            {
                n_act = n_act + (n_mid/n_norm)*integral_number(slope_2, n_origin_2, v_mid, v_int_2);
                v_act = v_act + (n_mid*v_mid/v_norm)*integral_volume(slope_2, n_origin_2, v_mid, v_int_2);
            }

            if (range[3])
            {
                n_act = n_act + (n_mid/n_norm)*integral_number(slope_2, n_origin_2, v_int_2, v_hi);
                v_act = v_act + (n_mid*v_mid/v_norm)*integral_volume(slope_2, n_origin_2, v_int_2, v_hi);
            }

            number_act[cb] = Kokkos::min(n_act, n_mid);
            for (int s=0; s<n_species; ++s)
                volume_act[cb][s] = Kokkos::min(v_act*(volume[ab][s]/(v_mid*n_mid)), volume[ab][s]);
        }

        // Move the activated particles from the aerosol to the cloud bins.
        for (int cb=cld_a_start; cb<cld_b_end; ++cb)
        {
            const int ab = parallel_aer(cb);

            number[ab] = Kokkos::max(0._tf, number[ab] - number_act[cb]);
            number[cb] = number[cb] + number_act[cb];

            for (int s=0; s<n_species; ++s)
            {
                volume[ab][s] = Kokkos::max(0._tf, volume[ab][s] - volume_act[cb][s]);
                volume[cb][s] = volume[cb][s] + volume_act[cb][s];
            }
        }
    }
}
