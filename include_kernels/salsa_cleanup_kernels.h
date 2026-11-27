/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), constrain_SALSA.f90, util.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Removal of small particles: evaporated cloud and rain back to aerosol, empty bins, and tiny aerosol.
 */

#pragma once

#include "salsa_types.h"
#include "salsa_constants.h"
#include "salsa_layout.h"


namespace Salsa
{
    // Fortran: util.f90: calc_correlation
    // Meant as the Pearson correlation of x and y, with two quirks kept: the cross term holds only the last
    // product, and r = 0 if a variance term is below epsilon (absolute).
    KOKKOS_INLINE_FUNCTION
    TF correlation(const Species_array& x, const Species_array& y, const int n)
    {
        TF sx = 0._tf;
        TF sy = 0._tf;
        TF sxx = 0._tf;
        TF syy = 0._tf;

        for (int s=0; s<n; ++s)
        {
            sx = sx + x[s];
            sy = sy + y[s];
        }
        for (int s=0; s<n; ++s)
            sxx = sxx + x[s]*x[s];
        for (int s=0; s<n; ++s)
            syy = syy + y[s]*y[s];

        const TF sxy = x[n-1]*y[n-1];
        const TF nf = static_cast<TF>(n);

        if (nf*sxx - sx*sx < eps || nf*syy - sy*sy < eps)
            return 0._tf;

        return (nf*sxy - sx*sy)/(Kokkos::sqrt(nf*sxx - sx*sx)*Kokkos::sqrt(nf*syy - sy*sy));
    }


    // Fortran: constrain_SALSA.f90: findDry4Wet
    // Aerosol bin for evaporating rain bin b: the 2a bin of its dry diameter, or the parallel bin when only that
    // one holds particles; if both do, the one whose dry composition correlates best.
    template<typename Number, typename Mass>
    KOKKOS_INLINE_FUNCTION
    int find_aerosol_bin_for_rain(
            const Number& number,
            const Mass& mass,
            const Array_1d<TF>& density,
            const Array_1d<TF>& diameter_low,
            const int b,
            const int n_species,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start)
    {
        const int n_dry = n_species - 1;

        TF d_dry = 0._tf;
        for (int s=0; s<n_dry; ++s)
            d_dry = d_dry + mass(b, s)/density(s);
        d_dry = d_dry/number(b);
        d_dry = Kokkos::pow(d_dry/pi6, 1._tf/3._tf);

        int n_larger = 0;
        for (int ba=aer_2a_start; ba<aer_2a_end; ++ba)
            if (d_dry > diameter_low(ba))
                ++n_larger;

        const int b_a = Kokkos::min(Kokkos::max(aer_2a_start + n_larger - 1, aer_2a_start), aer_2a_end-1);

        const int b_b = aer_2b_start + (b_a - aer_2a_start);

        if (number(b_b) <= nlim)
            return b_a;
        else if (number(b_a) <= nlim)
            return b_b;
        else
        {
            Species_array mass_rain;
            Species_array mass_a;
            Species_array mass_b;

            for (int s=0; s<n_dry; ++s)
            {
                mass_rain[s] = mass(b, s);
                mass_a[s] = mass(b_a, s);
                mass_b[s] = mass(b_b, s);
            }

            const TF r_a = correlation(mass_a, mass_rain, n_dry);
            const TF r_b = correlation(mass_b, mass_rain, n_dry);

            return (r_a < r_b) ? b_b : b_a;
        }
    }


    // Number and mass of bin b to zero (accessors return references to per-kg values).
    template<typename Number, typename Mass>
    KOKKOS_INLINE_FUNCTION
    void empty_bin(const Number& number, const Mass& mass, const int b, const int n_species)
    {
        number(b) = 0._tf;
        for (int s=0; s<n_species; ++s)
            mass(b, s) = 0._tf;
    }


    // Fortran: constrain_SALSA.f90: SALSA_diagnostics (part)
    // Empties bins with number but no mass; rain checks water only. Returns the water removed (kg kg-1).
    template<typename Number, typename Mass>
    KOKKOS_INLINE_FUNCTION
    TF empty_massless_bins(
            const Number& number,
            const Mass& mass,
            const int n_species,
            const int aer_1a_start,
            const int aer_2b_end,
            const int cld_a_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end)
    {
        const int iwa = n_species - 1;
        TF water_removed = 0._tf;

        auto mass_total = [&](const int b)
        {
            TF m = 0._tf;
            for (int s=0; s<n_species; ++s)
                m = m + mass(b, s);
            return m;
        };

        auto empty = [&](const int b)
        {
            water_removed += mass(b, iwa);
            empty_bin(number, mass, b, n_species);
        };

        for (int b=aer_1a_start; b<aer_2b_end; ++b)
            if (number(b) > 0._tf && mass_total(b) <= 0._tf)
                empty(b);

        for (int b=cld_a_start; b<cld_b_end; ++b)
            if (number(b) > 0._tf && mass_total(b) <= 0._tf)
                empty(b);

        for (int b=rain_start; b<rain_end; ++b)
            if (number(b) > 0._tf && mass(b, iwa) <= 0._tf)
                empty(b);

        return water_removed;
    }


    // Fortran: constrain_SALSA.f90: SALSA_diagnostics (part)
    // Removes aerosol below 0.1 nm dry diameter; SO4, OC, NO and NH to the gas phase (molecules per kg).
    // Returns the water removed (kg kg-1).
    template<typename Number, typename Mass, typename Gas_accessor>
    KOKKOS_INLINE_FUNCTION
    TF remove_tiny_aerosol(
            const Number& number,
            const Mass& mass,
            const Gas_accessor& gas,
            const Array_1d<TF>& density,
            const bool condense_gases,
            const int iso4,
            const int ioc,
            const int ino,
            const int inh,
            const int n_species,
            const int aer_1a_start,
            const int aer_2b_end)
    {
        const int n_dry = n_species - 1;
        const int iwa = n_species - 1;
        const TF volume_min = pi6*Kokkos::pow(1.e-10_tf, 3._tf);
        TF water_removed = 0._tf;

        for (int b=aer_1a_start; b<aer_2b_end; ++b)
        {
            if (number(b) > nlim)
            {
                TF volume = 0._tf;
                for (int s=0; s<n_dry; ++s)
                    volume = volume + mass(b, s)/density(s);
                volume = volume/number(b);

                if (volume < volume_min)
                {
                    auto to_gas = [&](const int s, const int g, const TF molar_mass_gas)
                    {
                        if (s >= 0 && condense_gases)
                            gas(g) = gas(g) + mass(b, s)/molar_mass_gas*avog;
                    };

                    to_gas(iso4, Gas::h2so4, msu);
                    to_gas(ioc, Gas::ocsv, moc);
                    to_gas(ino, Gas::hno3, mno);
                    to_gas(inh, Gas::nh3, mnh);

                    water_removed += mass(b, iwa);
                    empty_bin(number, mass, b, n_species);
                }
            }
        }

        return water_removed;
    }


    // Fortran: constrain_SALSA.f90: SALSA_diagnostics (part)
    // In subsaturated air, cloud and rain bins with small drops go to an aerosol bin, with their water
    // (accessors return references to per-kg values).
    template<typename Number, typename Mass>
    KOKKOS_INLINE_FUNCTION
    void evaporate_to_aerosol(
            const Number& number,
            const Mass& mass,
            const Array_1d<TF>& density,
            const Array_1d<TF>& molar_mass,
            const Array_1d<TF>& dissociation,
            const Array_1d<TF>& diameter_low,
            const TF T,
            const TF rh,
            const int n_species,
            const int aer_2a_start,
            const int aer_2a_end,
            const int aer_2b_start,
            const int cld_a_start,
            const int cld_a_end,
            const int cld_b_start,
            const int cld_b_end,
            const int rain_start,
            const int rain_end)
    {
        const int n_dry = n_species - 1;
        const int iwa = n_species - 1;
        const TF mass_threshold = 1.e-25_tf;           // Fortran: massTH

        // Moves number and mass of bin b to aerosol bin b_aer.
        auto move_to_aerosol = [&](const int b, const int b_aer)
        {
            number(b_aer) = number(b_aer) + number(b);
            number(b) = 0._tf;

            for (int s=0; s<n_species; ++s)
            {
                mass(b_aer, s) = mass(b_aer, s) + mass(b, s);
                mass(b, s) = 0._tf;
            }
        };

        // Wet diameter, and critical diameter from the Raoult and Kelvin terms.
        auto diameters = [&](const int b, TF& d_wet, TF& d_critical)
        {
            TF moles = 0._tf;
            for (int s=0; s<n_dry; ++s)
                moles = moles + dissociation(s)*mass(b, s)/molar_mass(s);
            moles = moles/number(b);

            const TF raoult = 6._tf*mwa*moles/(pi*rhowa);
            const TF kelvin = 4._tf*mwa*surfw0/(rg*rhowa*T);
            d_critical = Kokkos::sqrt(3._tf*raoult/kelvin);

            TF volume = 0._tf;
            for (int s=0; s<n_species; ++s)
                volume = volume + mass(b, s)/density(s);
            volume = volume/number(b);
            d_wet = Kokkos::pow(volume/pi6, 1._tf/3._tf);
        };

        // Evaporated cloud droplets to the parallel aerosol bin.
        for (int b=cld_a_start; b<cld_b_end; ++b)
        {
            if (number(b) > nlim && rh < 0.999_tf && mass(b, iwa) < 1.e-5_tf)
            {
                TF d_wet;
                TF d_critical;
                diameters(b, d_wet, d_critical);

                if (d_wet < Kokkos::min(0.7_tf*d_critical, 1.3e-5_tf) || mass(b, iwa) < mass_threshold*number(b))
                {
                    const int b_aer = (b < cld_a_end) ? aer_2a_start + (b - cld_a_start) : aer_2b_start + (b - cld_b_start);
                    move_to_aerosol(b, b_aer);
                }
            }
        }

        // Evaporated rain drops to the aerosol bin of their dry size and composition.
        for (int b=rain_start; b<rain_end; ++b)
        {
            if (number(b) > prlim && rh < 0.999_tf && mass(b, iwa) < 1.e-6_tf)
            {
                TF d_wet;
                TF d_critical;
                diameters(b, d_wet, d_critical);

                if (d_wet < Kokkos::min(0.2_tf*d_critical, 10.e-6_tf) || mass(b, iwa) < mass_threshold*number(b))
                {
                    const int b_aer = find_aerosol_bin_for_rain(
                            number,
                            mass,
                            density,
                            diameter_low,
                            b,
                            n_species,
                            aer_2a_start,
                            aer_2a_end,
                            aer_2b_start);
                    move_to_aerosol(b, b_aer);
                }
            }
        }
    }


    // Fortran: constrain_SALSA.f90: SALSA_diagnostics (without clipping, ice and charge)
    // On the fields of one grid point (per kg): bins with number but no mass are emptied; in subsaturated air,
    // cloud and rain bins with small drops go to an aerosol bin; aerosol below 0.1 nm dry diameter is removed,
    // with its volatile species to the gas phase.
    KOKKOS_INLINE_FUNCTION
    void remove_small_particles(
            const Array_4d<TF>& fields,
            const int k,
            const int j,
            const int i,
            const Array_1d<TF>& density,
            const Array_1d<TF>& molar_mass,
            const Array_1d<TF>& dissociation,
            const Array_1d<TF>& diameter_low,
            const TF T,
            const TF rh,
            const bool condense_gases,
            const int iso4,
            const int ioc,
            const int ino,
            const int inh,
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
            const int rain_end,
            const int n_bins)
    {
        auto number = [&](const int b) -> TF&
        {
            return fields(number_field(b, n_species), k, j, i);
        };

        auto mass = [&](const int b, const int s) -> TF&
        {
            return fields(mass_field(b, s, n_species), k, j, i);
        };

        auto gas = [&](const int g) -> TF&
        {
            return fields(gas_field(g, n_bins, n_species), k, j, i);
        };

        // Bins with number but no mass (their water is lost, as in the Fortran).
        empty_massless_bins(number, mass, n_species, aer_1a_start, aer_2b_end, cld_a_start, cld_b_end, rain_start, rain_end);

        evaporate_to_aerosol(
                number,
                mass,
                density,
                molar_mass,
                dissociation,
                diameter_low,
                T,
                rh,
                n_species,
                aer_2a_start,
                aer_2a_end,
                aer_2b_start,
                cld_a_start,
                cld_a_end,
                cld_b_start,
                cld_b_end,
                rain_start,
                rain_end);

        // Aerosol below 0.1 nm dry diameter (water lost, as in the Fortran).
        remove_tiny_aerosol(
                number, mass, gas, density, condense_gases, iso4, ioc, ino, inh, n_species, aer_1a_start, aer_2b_end);
    }
}
