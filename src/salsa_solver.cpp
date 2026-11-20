/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * The solver: one SALSA step for all grid points.
 */

#include "salsa_solver.h"


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
    }


    void Solver::exec(
            const Array_4d<TF>& fields,
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
        Kokkos::parallel_for("salsa_exec",
                Range_3d({kstart, jstart, istart}, {kend, jend, iend}),
                KOKKOS_LAMBDA(const int k, const int j, const int i)
                {
                    // The processes come here (step 5 on).
                });
    }
}
