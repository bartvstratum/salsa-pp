/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), mo_salsa_dynamics.f90, mo_salsa_coagulation_kernels.f90,
 * mo_salsa_coagulation_processes.f90, mo_particle_external_properties.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Coagulation: collisions within and between the particle categories.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"
#include "salsa_condensation_kernels.h"


namespace Salsa
{
    // Fortran: mo_particle_external_properties.f90: terminal_vel (aerosol, cloud and rain)
    KOKKOS_INLINE_FUNCTION
    TF terminal_velocity(
            const TF diam,                                  // m
            const TF rhop,                                  // Particle density (kg m-3)
            const TF rhoa,                                  // Air density (kg m-3)
            const TF visc,                                  // Dynamic viscosity of air (kg m-1 s-1)
            const TF beta)                                  // Cunningham correction factor
    {
        if (diam <= 19.e-6_tf)
        {
            // Small cloud droplets.
            return (diam*diam)*(rhop - rhoa)*grav*beta/(18._tf*visc);
        }
        else if (diam > 19.e-6_tf && diam < 1.07e-3_tf)
        {
            // Large cloud droplets and small rain drops.
            const TF n_davies = 4._tf/3._tf*Kokkos::pow(diam, 3._tf)*grav*(rhop - rhoa)*rhoa/(visc*visc);
            const TF x = Kokkos::log(n_davies);
            const TF x2 = x*x;
            const TF x3 = x2*x;
            const TF x4 = x2*x2;
            const TF x5 = x3*x2;
            const TF x6 = x3*x3;
            const TF y = -0.318657e1_tf + 0.992696_tf*x - 0.153193e-2_tf*x2 - 0.987059e-3_tf*x3
                       - 0.578878e-3_tf*x4 + 0.855176e-4_tf*x5 - 0.327815e-5_tf*x6;
            const TF reynolds = beta*Kokkos::exp(y);
            return visc/rhoa/diam*reynolds;
        }
        else
        {
            // Rain drops.
            const TF n_p = Kokkos::pow(surfw0, 3._tf)*(rhoa*rhoa)/(Kokkos::pow(visc, 4._tf)*(rhop - rhoa)*grav);
            const TF bond = 4._tf/3._tf*(rhop - rhoa)*grav/surfw0*(diam*diam);
            const TF x = Kokkos::log(bond*Kokkos::pow(n_p, 1._tf/6._tf));
            const TF x2 = x*x;
            const TF x3 = x2*x;
            const TF x4 = x2*x2;
            const TF x5 = x3*x2;
            const TF y = -0.500015e1_tf + 0.523778e1_tf*x - 0.204914e1_tf*x2 + 0.475294_tf*x3
                       - 0.542819e-1_tf*x4 + 0.238449e-2_tf*x5;
            const TF reynolds = Kokkos::pow(n_p, 1._tf/6._tf)*Kokkos::exp(y);
            return Kokkos::min(visc/rhoa/diam*reynolds, 9.2_tf);
        }
    }


    // Fortran: mo_salsa_coagulation_kernels.f90: coagc (kernel 2, without ice and charge)
    // Collision kernel (m3 s-1) of two particles with density rhowa: Brownian, convective enhancement,
    // gravitational collection and turbulence.
    KOKKOS_INLINE_FUNCTION
    TF coagulation_coefficient(
            const TF d_wet_1,
            const TF d_mid_1,
            const TF d_wet_2,
            const TF d_mid_2,
            const TF T,
            const TF p)
    {
        TF diam[2] = {d_wet_1, d_wet_2};
        TF mass[2] = {rhowa*pi6*(d_wet_1*d_wet_1*d_wet_1), rhowa*pi6*(d_wet_2*d_wet_2*d_wet_2)};

        // Bins of the same size (also the same bin): the first particle is taken 10 % larger.
        if (Kokkos::abs(d_mid_1 - d_mid_2)/d_mid_1 < 1.e-2_tf)
        {
            diam[0] = 1.1_tf*diam[1];
            mass[0] = 1.33_tf*mass[1];
        }

        const TF visc = (7.44523e-3_tf*Kokkos::sqrt(T*T*T))/(5093._tf*(T + 110.4_tf));
        const TF mfp = (1.656e-10_tf*T + 1.828e-8_tf)*pstand/p;    // Mean free path of air (m)

        TF beta[2];
        TF dfpart[2];
        TF mtvel[2];
        TF omega[2];

        for (int n=0; n<2; ++n)
        {
            const TF knud = 2._tf*mfp/diam[n];
            beta[n] = 1._tf + knud*(1.142_tf + 0.558_tf*Kokkos::exp(-0.999_tf/knud));
            dfpart[n] = beta[n]*boltz*T/(3._tf*pi*visc*diam[n]);
            mtvel[n] = Kokkos::sqrt((8._tf*boltz*T)/(pi*mass[n]));
            omega[n] = 8._tf*dfpart[n]/(pi*mtvel[n]);
        }

        const TF mdiam = 0.5_tf*(diam[0] + diam[1]);

        const TF flux_continuum = 4._tf*pi*mdiam*(dfpart[0] + dfpart[1]);
        const TF flux_free = pi*Kokkos::sqrt((mtvel[0]*mtvel[0]) + (mtvel[1]*mtvel[1]))*(mdiam*mdiam);

        TF tva[2];
        for (int n=0; n<2; ++n)
        {
            const TF sum = mdiam + omega[n];
            const TF sum_squares = mdiam*mdiam + omega[n]*omega[n];
            tva[n] = ((sum*sum*sum) - sum_squares*Kokkos::sqrt(sum_squares))/(3._tf*mdiam*omega[n]) - mdiam;
        }

        const TF fmdist = Kokkos::sqrt(tva[0]*tva[0] + tva[1]*tva[1]);     // Flux matching distance (m)

        const int lrg = (diam[0] >= diam[1]) ? 0 : 1;
        const int sml = 1 - lrg;

        const TF rhoa = p/(rd*T);
        const TF vkin = visc/rhoa;

        TF termv[2];
        TF reyn[2];
        TF schm[2];

        for (int n=0; n<2; ++n)
        {
            termv[n] = terminal_velocity(diam[n], rhowa, rhoa, visc, beta[n]);
            reyn[n] = diam[n]*termv[n]/vkin;
            schm[n] = vkin/dfpart[n];
        }

        const TF stok = 2._tf*termv[sml]*Kokkos::abs(termv[0] - termv[1])/(diam[lrg]*grav);

        const TF zbrown = flux_continuum/(mdiam/(mdiam + fmdist) + flux_continuum/flux_free);

        TF zbrconv = 0._tf;
        if (reyn[lrg] <= 1._tf)
            zbrconv = 0.45_tf*zbrown*(Kokkos::pow(reyn[lrg]*schm[sml], 1._tf/3._tf));
        else if (reyn[lrg] > 1._tf)
            zbrconv = 0.45_tf*zbrown*Kokkos::sqrt(reyn[lrg])*(Kokkos::pow(schm[sml], 1._tf/3._tf));

        const TF d_mean = 0.5_tf*(diam[0] + diam[1]);
        const TF ztshear = Kokkos::sqrt(8._tf*pi*eddy_dis/(15._tf*vkin))*(d_mean*d_mean*d_mean);
        const TF zturbinert = pi*Kokkos::pow(eddy_dis, 0.75_tf)/(grav*Kokkos::sqrt(Kokkos::sqrt(vkin)))
                            *(d_mean*d_mean)*Kokkos::abs(termv[0] - termv[1]);

        // Gravitational collection.
        const TF zea = (stok*stok)/((stok + 0.5_tf)*(stok + 0.5_tf));

        TF zev = 0._tf;
        if (stok > 1.214_tf)
        {
            zev = 0.75_tf*Kokkos::log(2._tf*stok)/(stok - 1.214_tf);
            zev = 1._tf/((1._tf + zev)*(1._tf + zev));
        }

        const TF zecoll = (60._tf*zev + zea*reyn[lrg])/(60._tf + reyn[lrg]);
        TF zgrav = zecoll*pi*(mdiam*mdiam);
        zgrav = zgrav*Kokkos::abs(termv[0] - termv[1]);

        return zbrown + zbrconv + Kokkos::sqrt(zgrav*zgrav + ztshear*ztshear + zturbinert*zturbinert);
    }


    // Fortran: mo_salsa_coagulation_kernels.f90: buildKernelSelf (one element)
    // Kernel between bins b1 and b2 of one category, from the state before coagulation: the bin with the lower
    // index first, zero if either bin holds less than its number limit.
    KOKKOS_INLINE_FUNCTION
    TF coagulation_kernel_self(
            const int b1,
            const int b2,
            const Bin_array& d_wet,
            const bool* has_particles,
            const Array_1d<TF>& diameter_mid,
            const TF T,
            const TF p)
    {
        const int lo = Kokkos::min(b1, b2);
        const int hi = Kokkos::max(b1, b2);

        if (!has_particles[lo] || !has_particles[hi])
            return 0._tf;

        return coagulation_coefficient(d_wet[lo], diameter_mid(lo), d_wet[hi], diameter_mid(hi), T, p);
    }


    // Fortran: mo_salsa_coagulation_processes.f90: accumulateSink (without ice and diagnostics)
    // Loss rate (s-1) of target bin b by collection by the bins [start, end).
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void accumulate_sink(
            TF& sink,
            const Kernel& kernel,
            const Bin_array& number,
            const int b,
            const int start,
            const int end,
            const TF multiplier)
    {
        for (int c=start; c<end; ++c)
            sink = sink + multiplier*kernel(b, c)*number[c];
    }


    // Fortran: mo_salsa_coagulation_processes.f90: accumulateSource (without diagnostics)
    // Volume gain rate per particle (m3 m-3 s-1 per m-3) of target bin b by collecting the bins [start, end).
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void accumulate_source(
            Species_array& source,
            const Kernel& kernel,
            const Bin_species_array& volume,
            const int b,
            const int start,
            const int end,
            const int n_species)
    {
        for (int c=start; c<end; ++c)
            for (int s=0; s<n_species; ++s)
                source[s] = source[s] + kernel(c, b)*volume[c][s];
    }


    // Fortran: mo_salsa_coagulation_processes.f90: accumulateSinkReverse
    // As accumulate_sink, with the smaller category as collector: loss of target bin b to the bins [start, end)
    // of a category that comes first in the kernel.
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void accumulate_sink_reverse(
            TF& sink,
            const Kernel& kernel,
            const Bin_array& number,
            const int b,
            const int start,
            const int end)
    {
        for (int c=start; c<end; ++c)
            sink = sink + kernel(c, b)*number[c];
    }


    // Fortran: mo_salsa_coagulation_processes.f90: accumulateSourceReverse
    // As accumulate_source, with the collected category second in the kernel.
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void accumulate_source_reverse(
            Species_array& source,
            const Kernel& kernel,
            const Bin_species_array& volume,
            const int b,
            const int start,
            const int end,
            const int n_species)
    {
        for (int c=start; c<end; ++c)
            for (int s=0; s<n_species; ++s)
                source[s] = source[s] + kernel(b, c)*volume[c][s];
    }


    // Fortran: mo_salsa_coagulation_processes.f90: accumulatePrecipFormation (without ice and diagnostics)
    // Autoconversion: collisions of cloud bin b with the cloud bins [start, end) that give a drop larger than the
    // smallest rain drop go directly to the rain bin of that size; the others are regular coagulation.
    // Self collection if [start, end) is bin b itself.
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void accumulate_precip_formation(
            Species_array& source,
            TF& sink,
            Species_array& volume_sink_self,                // Cloud volume lost by self collection into rain
            Rain_species_array& volume_rain,                // Volume source of each rain bin
            Rain_array& number_rain,                        // Number source of each rain bin
            const Kernel& kernel,
            const Bin_array& number,
            const Bin_species_array& volume,
            const Bin_array& d_wet,
            const Array_1d<TF>& diameter_low,
            const TF fix_coag,
            const int b,
            const int start,
            const int end,
            const int n_species,
            const int rain_start,
            const int rain_end)
    {
        const bool self_collection = (b == start && b == end-1);

        for (int c=start; c<end; ++c)
        {
            // Diameter after the collision.
            const TF d_new = Kokkos::pow(d_wet[b]*d_wet[b]*d_wet[b] + d_wet[c]*d_wet[c]*d_wet[c], 1._tf/3._tf);

            // Rain bin of that diameter, counted as in UCLALES-SALSA.
            int n_larger = 0;
            for (int r=rain_start; r<rain_end; ++r)
                if (d_new > diameter_low(r))
                    ++n_larger;

            if (d_new > diameter_low(rain_start))
            {
                const int r = n_larger - 1;

                if (self_collection)
                {
                    for (int s=0; s<n_species; ++s)
                        volume_rain[r][s] = volume_rain[r][s] + volume[c][s]*number[b]*kernel(c, b)*fix_coag;

                    number_rain[r] = number_rain[r] + (0.5_tf*kernel(c, b)*number[c]*number[b])*fix_coag;

                    for (int s=0; s<n_species; ++s)
                        volume_sink_self[s] = volume_sink_self[s] + volume[c][s]*number[b]*kernel(c, b)*fix_coag;

                    sink = sink + 0.5_tf*kernel(c, b)*number[b];
                }
                else
                {
                    for (int s=0; s<n_species; ++s)
                        volume_rain[r][s] = volume_rain[r][s]
                                          + (volume[c][s]*number[b] + volume[b][s]*number[c])*kernel(c, b)*fix_coag;

                    number_rain[r] = number_rain[r] + (kernel(c, b)*number[c]*number[b])*fix_coag;

                    // The target bin loses these particles too.
                    sink = sink + kernel(b, c)*number[c];
                }
            }
            else
            {
                if (self_collection)
                    sink = sink + 0.5_tf*kernel(b, c)*number[c];
                else
                    for (int s=0; s<n_species; ++s)
                        source[s] = source[s] + kernel(c, b)*volume[c][s];
            }
        }
    }


    KOKKOS_INLINE_FUNCTION
    void apply_coagulation(
            Bin_array& number,
            Bin_species_array& volume,
            const Species_array& source,
            const TF sink,
            const TF sink_self,
            const TF dt,
            const int b,
            const int n_species)
    {
        for (int s=0; s<n_species; ++s)
            volume[b][s] = (volume[b][s] + dt*source[s]*number[b])/(1._tf + dt*sink);

        number[b] = number[b]/(1._tf + dt*(sink + sink_self));
    }


    // Fortran: mo_salsa_coagulation_processes.f90: applyCoagPrecipFormation (without ice)
    KOKKOS_INLINE_FUNCTION
    void apply_precip_formation(
            Bin_array& number,
            Bin_species_array& volume,
            const Species_array& volume_sink_self,
            const Rain_species_array& volume_rain,
            const Rain_array& number_rain,
            const TF dt,
            const int b,
            const int n_species,
            const int rain_start,
            const int rain_end)
    {
        for (int s=0; s<n_species; ++s)
            volume[b][s] = volume[b][s] - dt*volume_sink_self[s];

        for (int r=rain_start; r<rain_end; ++r)
        {
            for (int s=0; s<n_species; ++s)
                volume[r][s] = volume[r][s] + dt*volume_rain[r-rain_start][s];

            number[r] = number[r] + dt*number_rain[r-rain_start];
        }
    }


    // Fortran: mo_salsa_coagulation_processes.f90: coag_aero (without rain and ice)
    // Semi-implicit, bin by bin with the updated state of the bins before: a bin loses particles to larger bins
    // and to itself, and gains the volume of the smaller bins it collects.
    template<typename Kernel_aa, typename Kernel_ca>
    KOKKOS_INLINE_FUNCTION
    void coagulate_aerosol(
            Bin_array& number,
            Bin_species_array& volume,
            const Kernel_aa& kernel_aa,
            const Kernel_ca& kernel_ca,                     // Aerosol bin first, cloud bin second
            const Array_1d<TF>& number_limit,
            const TF dt,
            const bool coagulate_aer,
            const bool cld_collects_aer,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_a_end,
            const int cld_b_start,
            const int cld_b_end)
    {
        const TF one = 1._tf;
        const TF half = 0.5_tf;

        Species_array source;

        // Volume gained from cloud bins smaller than the aerosol bin (none with the current bins).
        auto source_from_cloud = [&](const int cld_a_last, const int cld_b_last, const int b)
        {
            if (cld_a_last >= cld_a_start)
                accumulate_source_reverse(source, kernel_ca, volume, b, cld_a_start, cld_a_last+1, n_species);
            if (cld_b_last >= cld_b_start)
                accumulate_source_reverse(source, kernel_ca, volume, b, cld_b_start, cld_b_last+1, n_species);
        };

        // Regime 1a.
        for (int b=aer_1a_start; b<aer_1a_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            // Corresponding cloud bins (or the first ones) in regimes a and b.
            int b_cld_a = cld_a_start + Kokkos::max(b - aer_2a_start, 0);
            int b_cld_b = cld_b_start + (b_cld_a - cld_a_start);

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            if (coagulate_aer)
            {
                accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2b_end, one);
                accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);
            }

            // Collection by larger and equal cloud bins.
            if (cld_collects_aer)
            {
                accumulate_sink(sink, kernel_ca, number, b, b_cld_a, cld_a_end, one);
                accumulate_sink(sink, kernel_ca, number, b, b_cld_b, cld_b_end, one);
            }

            if (coagulate_aer && b > aer_1a_start)
                accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b, n_species);

            if (cld_collects_aer)
            {
                b_cld_a = cld_a_start + (b - aer_2a_start) - 1;
                b_cld_b = cld_b_start + (b_cld_a - cld_a_start);
                source_from_cloud(b_cld_a, b_cld_b, b);
            }

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }

        // Regime 2a: also collected by the larger or equal 2b bins, and collects the smaller or equal 2b bins.
        for (int b=aer_2a_start; b<aer_2a_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            const int b_2b = aer_2b_start + (b - aer_2a_start);
            int b_cld_a = cld_a_start + Kokkos::max(b - aer_2a_start, 0);
            int b_cld_b = cld_b_start + Kokkos::max(b - aer_2a_start, 0);

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            if (coagulate_aer)
            {
                if (b < aer_2a_end-1)
                {
                    accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2a_end, one);
                    accumulate_sink(sink, kernel_aa, number, b, b_2b, aer_2b_end, one);
                }
                accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);
            }

            if (cld_collects_aer)
            {
                accumulate_sink(sink, kernel_ca, number, b, b_cld_a, cld_a_end, one);
                accumulate_sink(sink, kernel_ca, number, b, b_cld_b, cld_b_end, one);
            }

            if (coagulate_aer)
            {
                accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b, n_species);
                accumulate_source(source, kernel_aa, volume, b, aer_2b_start, b_2b+1, n_species);
            }

            if (cld_collects_aer)
            {
                b_cld_a = cld_a_start + (b - aer_2a_start) - 1;
                b_cld_b = cld_b_start + (b_cld_a - cld_a_start);
                source_from_cloud(b_cld_a, b_cld_b, b);
            }

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }

        // Regime 2b, the same with the roles of 2a and 2b swapped; the 1a bins are collected as part of 2a.
        for (int b=aer_2b_start; b<aer_2b_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            const int b_2a = aer_2a_start + (b - aer_2b_start);
            int b_cld_a = cld_a_start + Kokkos::max(b - aer_2b_start, 0);
            int b_cld_b = cld_b_start + Kokkos::max(b - aer_2b_start, 0);

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            if (coagulate_aer)
            {
                if (b < aer_2b_end-1)
                {
                    accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2b_end, one);
                    accumulate_sink(sink, kernel_aa, number, b, b_2a, aer_2a_end, one);
                }
                accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);
            }

            if (cld_collects_aer)
            {
                accumulate_sink(sink, kernel_ca, number, b, b_cld_a, cld_a_end, one);
                accumulate_sink(sink, kernel_ca, number, b, b_cld_b, cld_b_end, one);
            }

            if (coagulate_aer)
            {
                accumulate_source(source, kernel_aa, volume, b, aer_2b_start, b, n_species);
                accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b_2a+1, n_species);
            }

            if (cld_collects_aer)
            {
                b_cld_b = cld_b_start + (b - aer_2b_start) - 1;
                b_cld_a = cld_a_start + (b_cld_b - cld_b_start);
                source_from_cloud(b_cld_a, b_cld_b, b);
            }

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }
    }


    // Fortran: mo_salsa_coagulation_processes.f90: coag_cloud (without rain and ice)
    // As coag_aero. Cloud bins of regime a and b are parallel to aerosol 2a and 2b; a cloud bin collects the
    // aerosol up to its parallel bins, and is collected by the larger aerosol bins (reverse collection).
    // With autoconversion (mode 1), cloud-cloud collisions that form drops of rain size go to the rain bins.
    template<typename Kernel_cc, typename Kernel_ca>
    KOKKOS_INLINE_FUNCTION
    void coagulate_cloud(
            Bin_array& number,
            Bin_species_array& volume,
            const Kernel_cc& kernel_cc,
            const Kernel_ca& kernel_ca,                     // Aerosol bin first, cloud bin second
            const Bin_array& d_wet,
            const Array_1d<TF>& diameter_low,
            const Array_1d<TF>& number_limit,
            const TF dt,
            const bool autoconversion,
            const bool coagulate_cld,
            const bool cld_collects_aer,
            const int n_species,
            const int aer_1a_start,
            const int aer_2a_start,
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
        const TF one = 1._tf;
        const TF half = 0.5_tf;
        const TF fix_coag_min = 0.1_tf;

        Species_array source;
        Species_array volume_sink_self;
        Rain_species_array volume_rain;
        Rain_array number_rain;

        // Correction of the precipitation formation for the collisions within one time step.
        auto fix_coag_self = [&](const int b)
        {
            return Kokkos::max(1._tf - dt*(0.5_tf*kernel_cc(b, b)*(number[b]*number[b])/(number[b])), fix_coag_min);
        };

        auto fix_coag_range = [&](const int b, const int start, const int end)
        {
            TF sum = 0._tf;
            for (int c=start; c<end; ++c)
                sum = sum + kernel_cc(b, c)*number[c];

            return Kokkos::max(1._tf - dt*sum, fix_coag_min);
        };

        auto precip_formation = [&](TF& sink, const int b, const int start, const int end, const TF fix_coag)
        {
            accumulate_precip_formation(
                    source, sink, volume_sink_self, volume_rain, number_rain, kernel_cc, number, volume, d_wet,
                    diameter_low, fix_coag, b, start, end, n_species, rain_start, rain_end);
        };

        // Sinks and sources of cloud bin b in regime [own_start, own_end), with the parallel regime
        // [other_start, other_end); i is the index within the regime.
        auto coagulate_bin = [&](
                const int b, const int i, const int own_start, const int own_end, const int other_start, const int other_end)
        {
            const int b_other = other_start + i;

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
            {
                source[s] = 0._tf;
                volume_sink_self[s] = 0._tf;
            }
            for (int r=0; r<rain_end-rain_start; ++r)
            {
                for (int s=0; s<n_species; ++s)
                    volume_rain[r][s] = 0._tf;
                number_rain[r] = 0._tf;
            }

            // Collection by larger aerosol bins in 2a and 2b.
            if (cld_collects_aer && b < own_end-1)
            {
                accumulate_sink_reverse(sink, kernel_ca, number, b, aer_2a_start + i + 1, aer_2a_end);
                accumulate_sink_reverse(sink, kernel_ca, number, b, aer_2b_start + i + 1, aer_2b_end);
            }

            // Collection by larger cloud bins in both regimes, and self collection.
            if (coagulate_cld)
            {
                if (b < own_end-1)
                {
                    accumulate_sink(sink, kernel_cc, number, b, b+1, own_end, one);
                    accumulate_sink(sink, kernel_cc, number, b, b_other, other_end, one);
                }

                if (autoconversion)
                    precip_formation(sink_self, b, b, b+1, fix_coag_self(b));
                else
                    accumulate_sink(sink_self, kernel_cc, number, b, b, b+1, half);
            }

            // Volume gained from smaller and equal aerosol bins.
            if (cld_collects_aer)
            {
                accumulate_source(source, kernel_ca, volume, b, aer_1a_start, aer_2a_start + i + 1, n_species);
                accumulate_source(source, kernel_ca, volume, b, aer_2b_start, aer_2b_start + i + 1, n_species);
            }

            // Volume gained from smaller cloud bins in the own regime, and smaller or equal in the other.
            if (coagulate_cld && b > own_start)
            {
                if (autoconversion)
                {
                    precip_formation(sink, b, own_start, b, fix_coag_range(b, own_start, own_end));
                    precip_formation(sink, b, other_start, b_other+1, fix_coag_range(b, other_start, other_end));
                }
                else
                {
                    accumulate_source(source, kernel_cc, volume, b, own_start, b, n_species);
                    accumulate_source(source, kernel_cc, volume, b, other_start, b_other+1, n_species);
                }
            }

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);

            if (autoconversion)
                apply_precip_formation(
                        number, volume, volume_sink_self, volume_rain, number_rain, dt, b, n_species, rain_start, rain_end);
        };

        // Regime a.
        for (int b=cld_a_start; b<cld_a_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            coagulate_bin(b, b - cld_a_start, cld_a_start, cld_a_end, cld_b_start, cld_b_end);
        }

        // Regime b, the same with the roles of a and b swapped.
        for (int b=cld_b_start; b<cld_b_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            coagulate_bin(b, b - cld_b_start, cld_b_start, cld_b_end, cld_a_start, cld_a_end);
        }
    }


    // Fortran: mo_salsa_dynamics.f90: coagulation (without rain coagulation and ice; kernel mode 1)
    KOKKOS_INLINE_FUNCTION
    void coagulate(
            Bin_array& number,
            Bin_species_array& volume,
            Bin_array& d_wet,
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& diameter_limit,
            const Array_1d<TF>& diameter_low,
            const Array_1d<TF>& number_limit,
            const TF T,
            const TF p,
            const TF dt,
            const bool autoconversion,                      // Mode 1, part of cloud coagulation
            const bool coagulate_aer,
            const bool coagulate_cld,
            const bool cld_collects_aer,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_a_end,
            const int cld_b_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end,
            const int n_bins)
    {
        update_wet_diameter(number, volume, d_wet, diameter_mid, diameter_limit, number_limit, n_species, n_bins);

        // The kernels use the state before coagulation; Fortran: update_coagulation_kernels
        bool has_particles[max_bins];
        for (int b=0; b<n_bins; ++b)
            has_particles[b] = !(number[b] < number_limit(b));

        // Within a category (aerosol or cloud).
        auto kernel_self = [&](const int b1, const int b2)
        {
            return coagulation_kernel_self(b1, b2, d_wet, has_particles, diameter_mid, T, p);
        };

        // Fortran: buildKernel (one element): aerosol bin b_aer and cloud bin b_cld.
        auto kernel_ca = [&](const int b_aer, const int b_cld)
        {
            if (!has_particles[b_aer] || !has_particles[b_cld])
                return 0._tf;

            return coagulation_coefficient(
                    d_wet[b_aer], diameter_mid(b_aer), d_wet[b_cld], diameter_mid(b_cld), T, p);
        };

        if (coagulate_aer || cld_collects_aer)
            coagulate_aerosol(
                    number,
                    volume,
                    kernel_self,
                    kernel_ca,
                    number_limit,
                    dt,
                    coagulate_aer,
                    cld_collects_aer,
                    n_species,
                    aer_1a_start,
                    aer_1a_end,
                    aer_2a_start,
                    aer_2a_end,
                    aer_2b_start,
                    aer_2b_end,
                    cld_a_start,
                    cld_a_end,
                    cld_b_start,
                    cld_b_end);

        if (coagulate_cld || cld_collects_aer)
            coagulate_cloud(
                    number,
                    volume,
                    kernel_self,
                    kernel_ca,
                    d_wet,
                    diameter_low,
                    number_limit,
                    dt,
                    autoconversion,
                    coagulate_cld,
                    cld_collects_aer,
                    n_species,
                    aer_1a_start,
                    aer_2a_start,
                    aer_2a_end,
                    aer_2b_start,
                    aer_2b_end,
                    cld_a_start,
                    cld_a_end,
                    cld_b_start,
                    cld_b_end,
                    rain_start,
                    rain_end);

        // Small negative values from rounding.
        for (int b=0; b<n_bins; ++b)
        {
            number[b] = Kokkos::max(0._tf, number[b]);
            for (int s=0; s<n_species; ++s)
                volume[b][s] = Kokkos::max(0._tf, volume[b][s]);
        }
    }
}
