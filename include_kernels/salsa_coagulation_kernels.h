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


    // Fortran: mo_salsa_coagulation_processes.f90: applyCoag
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


    // Fortran: mo_salsa_coagulation_processes.f90: coag_aero (aerosol-aerosol only)
    // Semi-implicit, bin by bin with the updated state of the bins before: a bin loses particles to larger bins
    // and to itself, and gains the volume of the smaller bins it collects.
    template<typename Kernel>
    KOKKOS_INLINE_FUNCTION
    void coagulate_aerosol(
            Bin_array& number,
            Bin_species_array& volume,
            const Kernel& kernel_aa,
            const Array_1d<TF>& number_limit,
            const TF dt,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end)
    {
        const TF one = 1._tf;
        const TF half = 0.5_tf;

        Species_array source;

        // Regime 1a.
        for (int b=aer_1a_start; b<aer_1a_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2b_end, one);
            accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);

            if (b > aer_1a_start)
                accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b, n_species);

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }

        // Regime 2a: also collected by the larger or equal 2b bins, and collects the smaller or equal 2b bins.
        for (int b=aer_2a_start; b<aer_2a_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            const int b_2b = aer_2b_start + (b - aer_2a_start);

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            if (b < aer_2a_end-1)
            {
                accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2a_end, one);
                accumulate_sink(sink, kernel_aa, number, b, b_2b, aer_2b_end, one);
            }
            accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);

            accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b, n_species);
            accumulate_source(source, kernel_aa, volume, b, aer_2b_start, b_2b+1, n_species);

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }

        // Regime 2b, the same with the roles of 2a and 2b swapped; the 1a bins are collected as part of 2a.
        for (int b=aer_2b_start; b<aer_2b_end; ++b)
        {
            if (number[b] < number_limit(b))
                continue;

            const int b_2a = aer_2a_start + (b - aer_2b_start);

            TF sink = 0._tf;
            TF sink_self = 0._tf;
            for (int s=0; s<n_species; ++s)
                source[s] = 0._tf;

            if (b < aer_2b_end-1)
            {
                accumulate_sink(sink, kernel_aa, number, b, b+1, aer_2b_end, one);
                accumulate_sink(sink, kernel_aa, number, b, b_2a, aer_2a_end, one);
            }
            accumulate_sink(sink_self, kernel_aa, number, b, b, b+1, half);

            accumulate_source(source, kernel_aa, volume, b, aer_2b_start, b, n_species);
            accumulate_source(source, kernel_aa, volume, b, aer_1a_start, b_2a+1, n_species);

            apply_coagulation(number, volume, source, sink, sink_self, dt, b, n_species);
        }
    }


    // Fortran: mo_salsa_dynamics.f90: coagulation (aerosol-aerosol only, kernel mode 1)
    KOKKOS_INLINE_FUNCTION
    void coagulate(
            Bin_array& number,
            Bin_species_array& volume,
            Bin_array& d_wet,
            const Array_1d<TF>& diameter_mid,
            const Array_1d<TF>& diameter_limit,
            const Array_1d<TF>& number_limit,
            const TF T,
            const TF p,
            const TF dt,
            const bool coagulate_aer,
            const int n_species,
            const int aer_1a_start,
            const int aer_1a_end,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int aer_2b_end,
            const int n_bins)
    {
        update_wet_diameter(number, volume, d_wet, diameter_mid, diameter_limit, number_limit, n_species, n_bins);

        // The kernels use the state before coagulation; Fortran: update_coagulation_kernels
        bool has_particles[max_bins];
        for (int b=0; b<n_bins; ++b)
            has_particles[b] = !(number[b] < number_limit(b));

        auto kernel_aa = [&](const int b1, const int b2)
        {
            return coagulation_kernel_self(b1, b2, d_wet, has_particles, diameter_mid, T, p);
        };

        if (coagulate_aer)
            coagulate_aerosol(
                    number,
                    volume,
                    kernel_aa,
                    number_limit,
                    dt,
                    n_species,
                    aer_1a_start,
                    aer_1a_end,
                    aer_2a_start,
                    aer_2a_end,
                    aer_2b_start,
                    aer_2b_end);

        // Small negative values from rounding.
        for (int b=0; b<n_bins; ++b)
        {
            number[b] = Kokkos::max(0._tf, number[b]);
            for (int s=0; s<n_species; ++s)
                volume[b][s] = Kokkos::max(0._tf, volume[b][s]);
        }
    }
}
