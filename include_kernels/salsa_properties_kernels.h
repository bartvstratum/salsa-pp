/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_properties.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Particle properties.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"


namespace Salsa
{
    // Fortran: mo_salsa_properties.f90: equilibration
    // Aerosol water in equilibrium with the relative humidity: every call for regime 1a (SO4, OC, NO), and for
    // regimes 2a and 2b (SO4, OC, NO and sea salt) only for initialisation. NO and NH add their volume.
    KOKKOS_INLINE_FUNCTION
    void equilibrate_water(
            const Bin_array& number,                        // m-3
            Bin_species_array& volume,                      // m3 m-3
            Bin_array& d_wet,                               // m
            Bin_array& core,                                // m3
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& density,
            const Array_1d<TF>& molar_mass,
            const TF rh,
            const TF T,
            const bool init,
            const int iso4,                                 // Species indices, -1 if not used
            const int ioc,
            const int ibc,
            const int idu,
            const int iss,
            const int ino,
            const int inh,
            const int iwa,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2b_end)
    {
        const int max_iterations = 100;

        auto volume_of = [&](const int b, const int s)
        {
            return (s >= 0) ? volume[b][s] : 0._tf;
        };

        auto ratio_of = [&](const int s)
        {
            return (s >= 0) ? density(s)/molar_mass(s) : 0._tf;
        };

        const TF ratio_so4 = ratio_of(iso4);
        const TF ratio_oc = ratio_of(ioc);
        const TF ratio_ss = ratio_of(iss);
        const TF ratio_no = ratio_of(ino);

        // Binary molalities (mol kg-1) of sulphate and organic carbon.
        auto molality_so4 = [](const TF aw)
        {
            const TF aw2 = aw*aw;
            const TF aw3 = aw2*aw;
            const TF aw4 = aw2*aw2;
            return 1.1065495e+2_tf - 3.6759197e+2_tf*aw + 5.0462934e+2_tf*aw2 - 3.1543839e+2_tf*aw3 + 6.770824e+1_tf*aw4;
        };

        auto molality_oc = [](const TF aw)
        {
            return 1._tf/(aw*mwa) - 1._tf/mwa;
        };

        // Binary molality of nitric acid.
        auto molality_no = [](const TF aw)
        {
            const TF aw2 = aw*aw;
            const TF aw3 = aw2*aw;
            const TF aw4 = aw2*aw2;
            const TF aw5 = aw3*aw2;
            const TF aw6 = aw3*aw3;
            const TF aw7 = aw4*aw3;
            return 2.306844303e+1_tf - 3.563608869e+1_tf*aw - 6.210577919e+1_tf*aw2 + 5.510176187e+2_tf*aw3
                 - 1.460055286e+3_tf*aw4 + 1.894467542e+3_tf*aw5 - 1.220611402e+3_tf*aw6 + 3.098597737e+2_tf*aw7;
        };

        // Regime 1a. The Kelvin term carries over from one bin to the next, as in the Fortran.
        TF ke = 1.001_tf;

        for (int b=aer_1a_start; b<aer_1a_end; ++b)
        {
            if (number[b] > nlim)
            {
                const TF v_so4 = volume_of(b, iso4);
                const TF v_oc = volume_of(b, ioc);
                const TF v_no = volume_of(b, ino);
                const TF core_b = (0._tf + v_so4/number[b]) + v_oc/number[b];
                const TF v_no_nh = v_no/number[b] + volume_of(b, inh)/number[b];    // NO and NH in one particle

                const TF rh_b = Kokkos::min(Kokkos::max(rh, 0.05_tf), 0.98_tf);

                TF d = d_wet[b];
                TF d_old = 1._tf;
                TF lwc = 0._tf;
                int count = 0;

                while (Kokkos::abs(d/d_old - 1._tf) > 1.e-2_tf)
                {
                    d_old = Kokkos::max(d, 1.e-20_tf);

                    const TF aw = rh_b/ke;

                    lwc = (v_so4*ratio_so4)/molality_so4(aw) + epsoc*v_oc*ratio_oc/molality_oc(aw)
                        + (v_no*ratio_no)/molality_no(aw);

                    d = Kokkos::pow(lwc/number[b]/rhowa/pi6 + v_no_nh/pi6 + core_b/pi6, 1._tf/3._tf);
                    ke = Kokkos::exp(2._tf*surfw0*mvwa/(boltz*T*d));

                    ++count;
                    if (count > max_iterations)
                        Kokkos::abort("equilibrate_water: no convergence in regime 1a");
                }

                volume[b][iwa] = lwc/rhowa;

                if (init)
                {
                    d_wet[b] = d;
                    core[b] = core_b;
                }
            }
            else if (init)
            {
                const TF dm = diameter_mid(b);
                d_wet[b] = dm;
                core[b] = pi6*(dm*dm*dm);
            }
        }

        if (!init)
            return;

        // Regimes 2a and 2b.
        for (int b=aer_2a_start; b<aer_2b_end; ++b)
        {
            if (number[b] > nlim)
            {
                const TF v_so4 = volume_of(b, iso4);
                const TF v_oc = volume_of(b, ioc);
                const TF v_ss = volume_of(b, iss);
                const TF v_no = volume_of(b, ino);
                const TF v_no_nh = v_no/number[b] + volume_of(b, inh)/number[b];    // NO and NH in one particle

                const TF core_b = ((((0._tf + v_so4/number[b]) + v_oc/number[b]) + volume_of(b, ibc)/number[b])
                        + volume_of(b, idu)/number[b]) + v_ss/number[b];

                const TF rh_b = Kokkos::min(Kokkos::max(rh, 0.37_tf), 0.98_tf);

                TF ke_b = 1.02_tf;
                TF d = d_wet[b];
                TF d_old = 1._tf;
                TF lwc = 0._tf;
                int count = 0;
                TF change = 0._tf;                          // Relative change of d in the previous iteration

                while (Kokkos::abs(d/d_old - 1._tf) > 1.e-12_tf)
                {
                    d_old = Kokkos::max(d, 1.e-20_tf);

                    const TF aw = rh_b/ke_b;
                    const TF aw2 = aw*aw;
                    const TF aw3 = aw2*aw;
                    const TF aw4 = aw2*aw2;

                    // Binary molality of sea salt (NaCl).
                    const TF molality_ss =
                            5.875248e+1_tf - 1.8781997e+2_tf*aw + 2.7211377e+2_tf*aw2 - 1.8458287e+2_tf*aw3 + 4.153689e+1_tf*aw4;

                    lwc = (v_so4*ratio_so4)/molality_so4(aw) + epsoc*(v_oc*ratio_oc)/molality_oc(aw)
                        + (v_no*ratio_no)/molality_no(aw) + (v_ss*ratio_ss)/molality_ss;

                    d = Kokkos::pow(lwc/number[b]/rhowa/pi6 + v_no_nh/pi6 + core_b/pi6, 1._tf/3._tf);
                    ke_b = Kokkos::exp(2._tf*surfw0*mvwa/(boltz*T*d));

                    // Also stop at the rounding noise, where the change no longer decreases (1e-12 is below it in
                    // float). The first change is from the initial guess, not an iteration step: compare from the third.
                    const TF change_new = Kokkos::abs(d/d_old - 1._tf);
                    if (count >= 2 && change_new >= change)
                        break;
                    change = change_new;

                    ++count;
                    if (count > max_iterations)
                        Kokkos::abort("equilibrate_water: no convergence in regime 2");
                }

                volume[b][iwa] = lwc/rhowa;
                d_wet[b] = d;
                core[b] = core_b;
            }
            else
            {
                const TF dm = diameter_mid(b);
                d_wet[b] = dm;
                core[b] = pi6*(dm*dm*dm);
            }
        }
    }
}
