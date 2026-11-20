/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * The solver: one SALSA step for all grid points.
 */

#include <stdexcept>
#include <vector>

#include "salsa_solver.h"
#include "salsa_constants.h"
#include "salsa_properties_kernels.h"
#include "salsa_condensation_kernels.h"
#include "salsa_activation_kernels.h"


namespace
{
    using namespace Salsa;

    // Relative humidity of the initialisation call; Fortran: mo_salsa_driver.f90: init_rh
    constexpr TF init_rh = 0.3_tf;


    int species_index(const std::vector<std::string>& species, const std::string& name)
    {
        for (size_t s=0; s<species.size(); ++s)
            if (species[s] == name)
                return static_cast<int>(s);

        return -1;
    }


    Array_1d<TF> to_device(const std::string& name, const std::vector<TF>& values)
    {
        Array_1d<TF> array(name, values.size());
        auto array_h = Kokkos::create_mirror_view(array);
        for (size_t n=0; n<values.size(); ++n)
            array_h(n) = values[n];
        Kokkos::deep_copy(array, array_h);
        return array;
    }
}


namespace Salsa
{
    Solver::Solver(
            const std::array<TF, 4>& regime_limits,
            const std::array<int, 2>& n_bins_regime,
            const int n_bins_rain,
            const TF rain_diameter_low,
            const TF rain_volume_ratio,
            const std::vector<std::string>& dry_species,
            const int coagulation_mode,
            const int autoconversion_mode,
            const int activation_mode,
            const bool coagulate_aer,
            const bool coagulate_cld,
            const bool cld_collects_aer,
            const bool rain_collects_cld,
            const bool rain_collects_aer,
            const bool coagulate_rain,
            const bool condense_gases,
            const bool condense_water_aer,
            const bool condense_water_cld,
            const TF rh_limit) :
        layout(regime_limits, n_bins_regime, n_bins_rain, rain_diameter_low, rain_volume_ratio, dry_species),
        coagulation_mode(coagulation_mode),
        autoconversion_mode(autoconversion_mode),
        activation_mode(activation_mode),
        coagulate_aer(coagulate_aer),
        coagulate_cld(coagulate_cld),
        cld_collects_aer(cld_collects_aer),
        rain_collects_cld(rain_collects_cld),
        rain_collects_aer(rain_collects_aer),
        coagulate_rain(coagulate_rain),
        condense_gases(condense_gases),
        condense_water_aer(condense_water_aer),
        condense_water_cld(condense_water_cld),
        rh_limit(rh_limit),
        iso4(species_index(layout.species, "SO4")),
        ioc(species_index(layout.species, "OC")),
        ibc(species_index(layout.species, "BC")),
        idu(species_index(layout.species, "DU")),
        iss(species_index(layout.species, "SS")),
        iwa(species_index(layout.species, "H2O"))
    {
        if (coagulation_mode != 1 || autoconversion_mode != 1 || activation_mode != 1)
            throw std::runtime_error("Only mode 1 is supported for coagulation, autoconversion and activation");

        density = to_device("density", layout.density);
        molar_mass = to_device("molar_mass", layout.molar_mass);
        diameter_mid = to_device("diameter_mid", layout.diameter_mid);
        number_limit = to_device("number_limit", layout.number_limit);
        diameter_limit = to_device("diameter_limit", layout.diameter_limit);
        dissociation = to_device("dissociation", layout.dissociation);
        volume_ratio_low = to_device("volume_ratio_low", layout.volume_ratio_low);
        volume_ratio_high = to_device("volume_ratio_high", layout.volume_ratio_high);

        soluble = Array_1d<bool>("soluble", layout.soluble.size());
        auto soluble_h = Kokkos::create_mirror_view(soluble);
        for (size_t s=0; s<layout.soluble.size(); ++s)
            soluble_h(s) = layout.soluble[s];
        Kokkos::deep_copy(soluble, soluble_h);
    }


    void Solver::exec(
            const Array_4d<TF>& fields,
            const Array_4d<TF>& diameter_wet,
            const Array_3d<TF>& rv,
            const Array_3d<const TF>& T,
            const Array_3d<const TF>& p,
            const Array_3d<const TF>& rs,
            const Array_3d<const TF>& rho,
            const TF dt,
            const bool initialize,
            const bool do_coagulation,
            const bool update_coagulation_kernels,
            const bool do_condensation,
            const bool do_activation,
            const bool do_autoconversion,
            const bool do_update_distribution,
            const bool free_rh,
            const int istart,
            const int iend,
            const int jstart,
            const int jend,
            const int kstart,
            const int kend)
    {
        const int n_bins = layout.n_bins;
        const int n_species = layout.n_species;
        const int n_dry = n_species - 1;

        const int aer_1a_start = layout.aer_1a_start;
        const int aer_1a_end = layout.aer_1a_end;
        const int aer_2a_start = layout.aer_2a_start;
        const int aer_2a_end = layout.aer_2a_end;
        const int aer_2b_start = layout.aer_2b_start;
        const int aer_2b_end = layout.aer_2b_end;
        const int cld_a_start = layout.cld_a_start;
        const int cld_a_end = layout.cld_a_end;
        const int cld_b_start = layout.cld_b_start;
        const int cld_b_end = layout.cld_b_end;
        const int rain_start = layout.rain_start;
        const int rain_end = layout.rain_end;

        const int iso4 = this->iso4;
        const int ioc = this->ioc;
        const int ibc = this->ibc;
        const int idu = this->idu;
        const int iss = this->iss;
        const int iwa = this->iwa;

        const TF rh_limit = this->rh_limit;
        const bool do_condense_gases = this->condense_gases;
        const bool condense_water_aer = this->condense_water_aer;
        const bool condense_water_cld = this->condense_water_cld;

        const Array_1d<TF> density = this->density;
        const Array_1d<TF> molar_mass = this->molar_mass;
        const Array_1d<TF> diameter_mid = this->diameter_mid;
        const Array_1d<TF> number_limit = this->number_limit;
        const Array_1d<TF> diameter_limit = this->diameter_limit;
        const Array_1d<TF> dissociation = this->dissociation;
        const Array_1d<TF> volume_ratio_low = this->volume_ratio_low;
        const Array_1d<TF> volume_ratio_high = this->volume_ratio_high;
        const Array_1d<bool> soluble = this->soluble;

        Kokkos::parallel_for("salsa_exec",
                Range_3d({kstart, jstart, istart}, {kend, jend, iend}),
                KOKKOS_LAMBDA(const int k, const int j, const int i)
                {
                    Bin_array number;                       // m-3; Fortran: numc
                    Bin_species_array volume;               // m3 m-3; Fortran: volc
                    Bin_array d_wet;                        // m; Fortran: dwet
                    Bin_array core;                         // Dry volume of one particle (m3)
                    TF gas[Gas::count];                     // m-3

                    // Fortran: mo_salsa_driver.f90: run_SALSA (part)
                    const TF T_ijk = T(k, j, i);
                    const TF p_ijk = p(k, j, i);
                    const TF rs_ijk = rs(k, j, i);
                    const TF rho_ijk = rho(k, j, i);
                    TF rv_ijk = free_rh ? rv(k, j, i) : Kokkos::min(rv(k, j, i), rs_ijk*rh_limit);

                    const TF rv_old = rv_ijk;

                    for (int b=0; b<n_bins; ++b)
                    {
                        for (int s=0; s<n_species; ++s)
                            volume[b][s] = fields(mass_field(b, s, n_species), k, j, i)*rho_ijk/density(s);

                        number[b] = fields(number_field(b, n_species), k, j, i)*rho_ijk;
                        d_wet[b] = diameter_wet(b, k, j, i);

                        if (number[b] > number_limit(b))
                        {
                            TF volume_dry = 0._tf;
                            for (int s=0; s<n_dry; ++s)
                                volume_dry += volume[b][s];

                            core[b] = volume_dry/number[b];
                        }
                        else
                        {
                            const TF dm = diameter_mid(b);
                            core[b] = pi6*(dm*dm*dm);
                        }
                    }

                    for (int g=0; g<Gas::count; ++g)
                        gas[g] = fields(gas_field(g, n_bins, n_species), k, j, i)*rho_ijk;

                    if (initialize)
                        equilibrate_water(
                                number,
                                volume,
                                d_wet,
                                core,
                                diameter_mid,
                                density,
                                molar_mass,
                                init_rh,
                                T_ijk,
                                initialize,
                                iso4,
                                ioc,
                                ibc,
                                idu,
                                iss,
                                iwa,
                                aer_1a_start,
                                aer_1a_end,
                                aer_2a_start,
                                aer_2b_end);

                    // Fortran: mo_salsa_dynamics.f90: condensation
                    if (do_condensation)
                    {
                        if (do_condense_gases)
                            Salsa::condense_gases(
                                    number,
                                    volume,
                                    d_wet,
                                    gas[Gas::h2so4],
                                    gas[Gas::ocnv],
                                    gas[Gas::ocsv],
                                    diameter_mid,
                                    diameter_limit,
                                    number_limit,
                                    T_ijk,
                                    p_ijk,
                                    dt,
                                    iso4,
                                    ioc,
                                    n_species,
                                    aer_1a_start,
                                    aer_2a_start,
                                    aer_2b_end,
                                    cld_a_start,
                                    cld_b_end,
                                    rain_start,
                                    rain_end);

                        if (condense_water_aer || condense_water_cld)
                            condense_water(
                                    number,
                                    volume,
                                    d_wet,
                                    core,
                                    rv_ijk,
                                    diameter_mid,
                                    diameter_limit,
                                    number_limit,
                                    density,
                                    molar_mass,
                                    dissociation,
                                    T_ijk,
                                    p_ijk,
                                    rs_ijk,
                                    dt,
                                    condense_water_aer,
                                    condense_water_cld,
                                    iso4,
                                    ioc,
                                    ibc,
                                    idu,
                                    iss,
                                    iwa,
                                    n_species,
                                    aer_1a_start,
                                    aer_1a_end,
                                    aer_2a_start,
                                    aer_2b_end,
                                    cld_a_start,
                                    cld_b_end,
                                    rain_start,
                                    rain_end);
                    }

                    // Fortran: mo_salsa_cloud.f90: cloud_activation
                    if (do_activation)
                        activate(
                                number,
                                volume,
                                diameter_mid,
                                number_limit,
                                volume_ratio_low,
                                volume_ratio_high,
                                soluble,
                                T_ijk,
                                rv_ijk,
                                rs_ijk,
                                iwa,
                                n_species,
                                aer_1a_start,
                                aer_2a_start,
                                aer_2a_end,
                                aer_2b_start,
                                aer_2b_end,
                                cld_a_start,
                                cld_a_end,
                                cld_b_start,
                                cld_b_end);

                    // The other processes come here (step 6d on).

                    // Fortran: mo_salsa_driver.f90: run_SALSA (part)
                    // Fortran: constrain_SALSA.f90: tend_constrain2
                    // The change as a tendency, applied as in UCLALES-SALSA: a step lowers a field to at most 1e-8 of its
                    // old value, and the initialisation clips at zero instead.
                    auto update = [&](const int n, const TF value)
                    {
                        const TF old = fields(n, k, j, i);
                        const TF tendency = (value - old)/dt;

                        if (initialize)
                            fields(n, k, j, i) = Kokkos::max(old + dt*tendency, 0._tf);
                        else
                            fields(n, k, j, i) = old + dt*Kokkos::max((-1._tf + 1.e-8_tf)*old/dt, tendency);
                    };

                    for (int b=0; b<n_bins; ++b)
                    {
                        for (int s=0; s<n_species; ++s)
                            update(mass_field(b, s, n_species), volume[b][s]*density(s)/rho_ijk);

                        update(number_field(b, n_species), number[b]/rho_ijk);
                        diameter_wet(b, k, j, i) = d_wet[b];
                    }

                    if (do_condense_gases)
                        for (int g=0; g<Gas::count; ++g)
                            update(gas_field(g, n_bins, n_species), gas[g]/rho_ijk);

                    rv(k, j, i) += dt*((rv_ijk - rv_old)/dt);
                });
    }
}
