/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_dynamics.f90, classSection.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Condensation of gases and water vapour.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"
#include "salsa_properties_kernels.h"


namespace Salsa
{
    // Mass accommodation coefficient, 1 for all bins; Fortran: massacc
    constexpr TF massacc = 1._tf;


    // Fortran: classSection.f90: updateDiameter (wet diameter, limited)
    KOKKOS_INLINE_FUNCTION
    void update_wet_diameter(
            const Bin_array& number,
            const Bin_species_array& volume,
            Bin_array& d_wet,
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& diameter_limit,
            const Array_1d<TF>& number_limit,
            const int n_species,
            const int n_bins)
    {
        for (int b=0; b<n_bins; ++b)
        {
            if (number[b] > number_limit(b))
            {
                TF volume_total = 0._tf;
                for (int s=0; s<n_species; ++s)
                    volume_total += volume[b][s];

                d_wet[b] = Kokkos::pow(volume_total/number[b]/pi6, 1._tf/3._tf);
            }
            else
                d_wet[b] = diameter_mid(b);

            d_wet[b] = Kokkos::min(d_wet[b], diameter_limit(b));
        }
    }


    // Fortran: mo_salsa_dynamics.f90: acth2o
    // Water activity of a bin with the given water (mol m-3), assuming an activity coefficient of 1.
    KOKKOS_INLINE_FUNCTION
    TF water_activity(
            const Species_array& volume,                    // Volume per species of the bin (m3 m-3)
            const TF water,                                 // mol m-3
            const Array_1d<TF>& density,
            const Array_1d<TF>& molar_mass,
            const Array_1d<TF>& dissociation,
            const int n_dry)
    {
        TF solute = 0._tf;
        for (int s=0; s<n_dry; ++s)
            solute += dissociation(s)*volume[s]*density(s)/molar_mass(s);

        // A copy of eps: device code cannot bind a namespace constant to the reference argument of Kokkos::max.
        // A floor at the float epsilon (1.2e-7 mol m-3) would be as large as the water of small droplets.
        return Kokkos::max(0.1_tf, water/Kokkos::max(TF(eps), water + solute));
    }


    // Fortran: mo_salsa_dynamics.f90: condgas
    // Condensation of H2SO4 and organic vapours, without nucleation.
    KOKKOS_INLINE_FUNCTION
    void condense_gases(
            const Bin_array& number,                        // m-3
            Bin_species_array& volume,                      // m3 m-3
            Bin_array& d_wet,                               // m
            TF& h2so4,                                      // m-3
            TF& ocnv,                                       // Non-volatile organic vapour (m-3)
            TF& ocsv,                                       // Semi-volatile organic vapour (m-3)
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& diameter_limit,
            const Array_1d<TF>& number_limit,
            const TF T,
            const TF p,
            const TF dt,
            const int iso4,                                 // Species indices, -1 if not used
            const int ioc,
            const int n_species,
            const int aer_1a_start,
            const int aer_2a_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end)
    {
        update_wet_diameter(number, volume, d_wet, diameter_mid, diameter_limit, number_limit, n_species, rain_end);

        Bin_array dvol_oc;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            dvol_oc[b] = 0._tf;

        // Properties of air and the condensing gases.
        const TF visc = (7.44523e-3_tf*Kokkos::sqrt(T*T*T))/(5093._tf*(T + 110.4_tf));
        const TF dfvap = 5.1111e-10_tf*Kokkos::pow(T, 1.75_tf)*pstand/p;
        const TF mfp = 3._tf*dfvap*Kokkos::sqrt(pi*msu/(8._tf*rg*T));

        // Collision rate of molecules with particles (s-1), with the transition regime correction of Fuchs and
        // Sutugin (1971). The size of the molecule and particle diffusion only for the two smallest bins.
        Bin_array colrate;

        for (int b=aer_1a_start; b<aer_2b_end; ++b)
        {
            const bool smallest = (b < aer_1a_start + 2);

            const TF knud = smallest ? 2._tf*mfp/(d_wet[b] + d_sa) : 2._tf*mfp/d_wet[b];
            const TF beta = (knud + 1._tf)/(0.377_tf*knud + 1._tf + 4._tf/(3._tf*massacc)*(knud + knud*knud));

            if (number[b] > nlim)
            {
                if (smallest)
                {
                    const TF dfpart = boltz*T*beta/(3._tf*pi*visc*d_wet[b]);
                    colrate[b] = 2._tf*pi*(d_wet[b] + d_sa)*(dfvap + dfpart)*beta*number[b];
                }
                else
                    colrate[b] = 2._tf*pi*d_wet[b]*dfvap*beta*number[b];
            }
            else
                colrate[b] = 0._tf;
        }

        for (int b=cld_a_start; b<rain_end; ++b)
        {
            const TF knud = 2._tf*mfp/d_wet[b];
            TF beta = 1._tf + knud*(1.33_tf + (0.71_tf/knud))/(1._tf + (1._tf/knud));
            beta = 1._tf/beta;

            const TF limit = (b < rain_start) ? nlim : prlim;
            colrate[b] = (number[b] > limit) ? 2._tf*pi*d_wet[b]*dfvap*beta*number[b] : 0._tf;
        }

        // Condensation sinks (s-1).
        TF sum_aer = 0._tf;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            sum_aer += colrate[b];

        TF sum_aer_2 = 0._tf;
        for (int b=aer_2a_start; b<aer_2b_end; ++b)
            sum_aer_2 += colrate[b];

        TF sum_cld = 0._tf;
        for (int b=cld_a_start; b<cld_b_end; ++b)
            sum_cld += colrate[b];

        TF sum_rain = 0._tf;
        for (int b=rain_start; b<rain_end; ++b)
            sum_rain += colrate[b];

        const TF cs_tot = sum_aer + sum_cld + sum_rain;

        // Non-volatile organic vapour: all bins.
        if (ocnv > 1.e-10_tf && cs_tot > 1.e-30_tf && ioc >= 0)
        {
            const TF cvap_new = ocnv/(1._tf + dt*cs_tot);
            const TF dvap = ocnv - cvap_new;
            ocnv = cvap_new;

            for (int b=aer_1a_start; b<aer_2b_end; ++b)
            {
                dvol_oc[b] = colrate[b]/cs_tot*mvoc*dvap;
                volume[b][ioc] = volume[b][ioc] + dvol_oc[b];
            }

            for (int b=cld_a_start; b<rain_end; ++b)
                volume[b][ioc] = volume[b][ioc] + colrate[b]/cs_tot*mvoc*dvap;
        }

        // Semi-volatile organic vapour: regime 2, cloud and rain. As in the Fortran, the aerosol gets the non-volatile
        // change of this call again.
        const TF cs_ocsv = sum_aer_2 + sum_cld + sum_rain;

        if (ocsv > 1.e-10_tf && cs_ocsv > 1.e-30_tf && ioc >= 0)
        {
            const TF cvap_new = ocsv/(1._tf + dt*cs_ocsv);
            const TF dvap = ocsv - cvap_new;
            ocsv = cvap_new;

            for (int b=aer_2a_start; b<aer_2b_end; ++b)
                dvol_oc[b] = dvol_oc[b] + colrate[b]/cs_ocsv*mvoc*dvap;

            for (int b=aer_1a_start; b<aer_2b_end; ++b)
                volume[b][ioc] = volume[b][ioc] + dvol_oc[b];

            for (int b=cld_a_start; b<rain_end; ++b)
                volume[b][ioc] = volume[b][ioc] + colrate[b]/cs_ocsv*mvoc*dvap;
        }

        // Sulphuric acid: all bins.
        if (h2so4 > 1.e-10_tf && cs_tot > 1.e-30_tf && iso4 >= 0)
        {
            const TF cvap_new = h2so4/(1._tf + dt*cs_tot);
            const TF dvap = h2so4 - cvap_new;
            h2so4 = cvap_new;

            for (int b=aer_1a_start; b<aer_2b_end; ++b)
                volume[b][iso4] = volume[b][iso4] + colrate[b]/cs_tot*mvsu*dvap;

            for (int b=cld_a_start; b<rain_end; ++b)
                volume[b][iso4] = volume[b][iso4] + colrate[b]/cs_tot*mvsu*dvap;
        }
    }


    // Fortran: mo_salsa_dynamics.f90: gpparth2o
    // Condensation of water vapour on aerosol, cloud and rain, with the analytical predictor of Jacobson (2005) in
    // fixed substeps of 0.02 s. Regime 1a and, below RH 98 %, all aerosol are in equilibrium instead.
    KOKKOS_INLINE_FUNCTION
    void condense_water(
            const Bin_array& number,                        // m-3
            Bin_species_array& volume,                      // m3 m-3
            Bin_array& d_wet,                               // m
            Bin_array& core,                                // m3
            TF& rv,                                         // kg kg-1
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& diameter_limit,
            const Array_1d<TF>& number_limit,
            const Array_1d<TF>& density,
            const Array_1d<TF>& molar_mass,
            const Array_1d<TF>& dissociation,
            const TF T,
            const TF p,
            const TF rs,
            const TF dt,
            const bool condense_water_aer,
            const bool condense_water_cld,
            const int iso4,                                 // Species indices, -1 if not used
            const int ioc,
            const int ibc,
            const int idu,
            const int iss,
            const int iwa,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end)
    {
        const int n_dry = n_species - 1;
        const TF rh = rv/rs;

        // Equilibrium water of regime 1a, and of all aerosol below RH 98 %; the change goes to or from the vapour.
        TF lwc_old = 0._tf;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            lwc_old += volume[b][iwa];
        lwc_old = lwc_old*rhowa;

        const bool init_false = false;
        const bool init_true = true;

        equilibrate_water(
                number, volume, d_wet, core, diameter_mid, density, molar_mass, rh, T, init_false,
                iso4, ioc, ibc, idu, iss, iwa, aer_1a_start, aer_1a_end, aer_2a_start, aer_2b_end);

        if (rh < 0.98_tf || !condense_water_aer)
            equilibrate_water(
                    number, volume, d_wet, core, diameter_mid, density, molar_mass, rh, T, init_true,
                    iso4, ioc, ibc, idu, iss, iwa, aer_1a_start, aer_1a_end, aer_2a_start, aer_2b_end);

        TF lwc_new = 0._tf;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            lwc_new += volume[b][iwa];
        lwc_new = lwc_new*rhowa;

        rv = rv - (lwc_new - lwc_old)/(p*mair/(rg*T));

        // Nothing to do without particles that take up water.
        bool any_aer = false;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            if (number[b] > nlim)
                any_aer = true;

        bool any_cld = false;
        for (int b=cld_a_start; b<cld_b_end; ++b)
            if (number[b] > nlim)
                any_cld = true;

        bool any_rain = false;
        for (int b=rain_start; b<rain_end; ++b)
            if (number[b] > prlim)
                any_rain = true;

        if (!( (any_cld || (any_rain && condense_water_cld)) || (any_aer && rh > 0.98_tf && condense_water_aer) ))
            return;

        const TF rhoair = mair*p/(rg*T);

        TF dfh2o = (5._tf/(16._tf*avog*rhoair*1.e-3_tf*(3.11e-8_tf*3.11e-8_tf)))
                 * Kokkos::sqrt(rg*1.e7_tf*T*mair*1.e3_tf*(mwa + mair)*1.e3_tf/(2._tf*pi*mwa*1.e3_tf));
        dfh2o = dfh2o*1.e-4_tf;

        const TF mfph2o = 3._tf*dfh2o*Kokkos::sqrt(pi*mwa/(8._tf*rg*T));
        const TF thcond = 0.023807_tf + 7.1128e-5_tf*(T - 273.16_tf);

        update_wet_diameter(number, volume, d_wet, diameter_mid, diameter_limit, number_limit, n_species, rain_end);

        // Mass transfer coefficients, saturation ratios and Kelvin terms per bin.
        Bin_array kelvin;
        Bin_array cw_surf;
        Bin_array w_sat;
        Bin_array mt;

        for (int b=0; b<rain_end; ++b)
        {
            kelvin[b] = 1._tf;
            cw_surf[b] = 0._tf;
            w_sat[b] = 0._tf;
            mt[b] = 0._tf;

            const bool is_aer = (b < aer_2b_end);
            const bool active = is_aer ?
                    (number[b] > number_limit(b) && rh > 0.98_tf && condense_water_aer) :
                    (number[b] > number_limit(b) && condense_water_cld);

            if (!active)
                continue;

            const TF dwet = d_wet[b];
            const TF water = volume[b][iwa]*density(iwa)/molar_mass(iwa);
            const TF act = water_activity(volume[b], water, density, molar_mass, dissociation, n_dry);

            kelvin[b] = Kokkos::exp(4._tf*surfw0*mwa/(rg*T*rhowa*dwet));
            cw_surf[b] = rs*rhoair/mwa;
            w_sat[b] = act*kelvin[b];

            const TF knud = 2._tf*mfph2o/dwet;
            const TF beta = is_aer ?
                    (knud + 1._tf)/(0.377_tf*knud + 1._tf + 4._tf/(3._tf*massacc)*(knud + knud*knud)) :
                    (knud + 1._tf)/(0.377_tf*knud + 1._tf + 4._tf/(3._tf)*(knud + knud*knud));

            const TF hlp1 = number[b]*2._tf*pi*dwet*dfh2o*beta;
            const TF hlp2 = mwa*dfh2o*alv*w_sat[b]*cw_surf[b]/(thcond*T);
            const TF hlp3 = ((alv*mwa)/(rg*T)) - 1._tf;

            mt[b] = hlp1/(hlp2*hlp3 + 1._tf);
        }

        // Water in the vapour and in each bin (mol m-3).
        TF cw_vapour = rv*rhoair/mwa;
        Bin_array cw;
        for (int b=0; b<rain_end; ++b)
            cw[b] = volume[b][iwa]*rhowa/mwa;

        TF sum_aer = 0._tf;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            sum_aer += cw[b];

        TF sum_cld = 0._tf;
        for (int b=cld_a_start; b<cld_b_end; ++b)
            sum_cld += cw[b];

        TF sum_rain = 0._tf;
        for (int b=rain_start; b<rain_end; ++b)
            sum_rain += cw[b];

        const TF cw_total = cw_vapour + sum_aer + sum_cld + sum_rain;

        bool any_aer_section = false;
        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            if (number[b] > number_limit(b))
                any_aer_section = true;

        bool any_cld_section = false;
        for (int b=cld_a_start; b<cld_b_end; ++b)
            if (number[b] > number_limit(b))
                any_cld_section = true;

        bool any_rain_section = false;
        for (int b=rain_start; b<rain_end; ++b)
            if (number[b] > number_limit(b))
                any_rain_section = true;

        // Substeps of about 0.02 s that add up to dt; regime 1a is in equilibrium and does not take part.
        const TF adt_target = 2.e-2_tf;
        const int n_substeps = Kokkos::max(1, static_cast<int>(Kokkos::round(dt/adt_target)));
        const TF adt = dt/static_cast<TF>(n_substeps);
        const TF max_change = 1.e-1_tf;

        Bin_array cw_new;
        for (int b=0; b<rain_end; ++b)
            cw_new[b] = cw[b];

        TF cw_vapour_new = 0._tf;

        for (int n=0; n<n_substeps; ++n)
        {
            TF source_aer = 0._tf;
            TF sink_aer = 0._tf;
            for (int b=aer_2a_start; b<aer_2b_end; ++b)
            {
                source_aer += mt[b]*w_sat[b]*cw_surf[b];
                sink_aer += mt[b];
            }

            TF source_cld = 0._tf;
            TF sink_cld = 0._tf;
            for (int b=cld_a_start; b<cld_b_end; ++b)
            {
                source_cld += mt[b]*w_sat[b]*cw_surf[b];
                sink_cld += mt[b];
            }

            TF source_rain = 0._tf;
            TF sink_rain = 0._tf;
            for (int b=rain_start; b<rain_end; ++b)
            {
                source_rain += mt[b]*w_sat[b]*cw_surf[b];
                sink_rain += mt[b];
            }

            const TF hlp1 = cw_vapour + adt*(source_aer + source_cld + source_rain);
            const TF hlp2 = 1._tf + adt*(sink_aer + sink_cld + sink_rain);

            cw_vapour_new = Kokkos::min(hlp1/hlp2, cw_total);

            auto update = [&](const int b)
            {
                const TF change = adt*mt[b]*(cw_vapour_new - w_sat[b]*cw_surf[b]);
                cw_new[b] = cw[b] + Kokkos::min(Kokkos::max(change, -max_change*cw[b]), max_change*cw[b]);
                w_sat[b] = water_activity(volume[b], cw_new[b], density, molar_mass, dissociation, n_dry)*kelvin[b];
            };

            if (any_aer_section && rh > 0.98_tf)
                for (int b=aer_2a_start; b<aer_2b_end; ++b)
                    update(b);

            if (any_cld_section)
                for (int b=cld_a_start; b<cld_b_end; ++b)
                    update(b);

            if (any_rain_section)
                for (int b=rain_start; b<rain_end; ++b)
                    update(b);

            for (int b=aer_2a_start; b<rain_end; ++b)
                cw_new[b] = Kokkos::max(cw_new[b], 0._tf);

            TF sum_aer_new = 0._tf;
            for (int b=aer_1a_start; b<aer_2b_end; ++b)
                sum_aer_new += cw_new[b];

            TF sum_cld_new = 0._tf;
            for (int b=cld_a_start; b<cld_b_end; ++b)
                sum_cld_new += cw_new[b];

            TF sum_rain_new = 0._tf;
            for (int b=rain_start; b<rain_end; ++b)
                sum_rain_new += cw_new[b];

            cw_vapour_new = cw_total - sum_aer_new - sum_cld_new - sum_rain_new;

            for (int b=0; b<rain_end; ++b)
                cw[b] = cw_new[b];
            cw_vapour = cw_vapour_new;
        }

        rv = cw_vapour_new*mwa/rhoair;

        for (int b=0; b<rain_end; ++b)
            volume[b][iwa] = Kokkos::max(0._tf, cw_new[b]*mwa/rhowa);
    }
}
