/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * The solver: one SALSA step for all grid points.
 */

#include <vector>

#include "salsa_solver.h"


namespace
{
    using namespace Salsa;

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
        rh_limit(rh_limit)
    {
        density = to_device("density", layout.density);
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
        const int n_gases = Layout::n_gases;

        const TF rh_limit = this->rh_limit;
        const bool condense_gases = this->condense_gases;

        const Array_1d<TF> density = this->density;

        Kokkos::parallel_for("salsa_exec",
                Range_3d({kstart, jstart, istart}, {kend, jend, iend}),
                KOKKOS_LAMBDA(const int k, const int j, const int i)
                {
                    TF number[max_bins];                    // m-3; Fortran: numc
                    TF volume[max_bins][max_species];       // m3 m-3; Fortran: volc
                    TF d_wet[max_bins];                     // m; Fortran: dwet
                    TF gas[Layout::n_gases];                // m-3

                    // Fortran: mo_salsa_driver.f90: run_SALSA (part)
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
                    }

                    for (int g=0; g<n_gases; ++g)
                        gas[g] = fields(gas_field(g, n_bins, n_species), k, j, i)*rho_ijk;

                    // The processes come here (step 6b on).

                    // Fortran: mo_salsa_driver.f90: run_SALSA (part)
                    for (int b=0; b<n_bins; ++b)
                    {
                        for (int s=0; s<n_species; ++s)
                            fields(mass_field(b, s, n_species), k, j, i) = volume[b][s]*density(s)/rho_ijk;

                        fields(number_field(b, n_species), k, j, i) = number[b]/rho_ijk;
                        diameter_wet(b, k, j, i) = d_wet[b];
                    }

                    if (condense_gases)
                        for (int g=0; g<n_gases; ++g)
                            fields(gas_field(g, n_bins, n_species), k, j, i) = gas[g]/rho_ijk;

                    rv(k, j, i) += rv_ijk - rv_old;
                });
    }
}
