/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Test program: creates a solver with the bins and species of reference case C and calls it once on a small grid
 * with ghost cells, as MicroHH's fields.
 */

#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "salsa_solver.h"


int main(int argc, char* argv[])
{
    using namespace Salsa;

    Kokkos::ScopeGuard guard(argc, argv);

    // Bins and species of reference case C (all dry species, regimes 2a and 2b), master defaults otherwise.
    const std::array<TF, 4> regime_limits = {3.e-9_tf, 5.e-8_tf, 7.e-7_tf, 1.e-5_tf};
    const std::array<int, 2> n_bins_regime = {3, 7};
    const int n_bins_rain = 7;
    const TF rain_diameter_low = 50.e-6_tf;
    const TF rain_volume_ratio = 4.859562_tf;
    const std::vector<std::string> dry_species = {"SO4", "OC", "BC", "DU", "SS", "NO", "NH"};

    const int coagulation_mode = 1;
    const int autoconversion_mode = 1;
    const int activation_mode = 1;

    const bool coagulate_aer = true;
    const bool coagulate_cld = true;
    const bool cld_collects_aer = true;
    const bool rain_collects_cld = true;
    const bool rain_collects_aer = true;
    const bool coagulate_rain = true;

    const bool condense_gases = true;
    const bool condense_water_aer = true;
    const bool condense_water_cld = true;

    const TF rh_limit = 1.20_tf;

    Solver solver(
            regime_limits,
            n_bins_regime,
            n_bins_rain,
            rain_diameter_low,
            rain_volume_ratio,
            dry_species,
            coagulation_mode,
            autoconversion_mode,
            activation_mode,
            coagulate_aer,
            coagulate_cld,
            cld_collects_aer,
            rain_collects_cld,
            rain_collects_aer,
            coagulate_rain,
            condense_gases,
            condense_water_aer,
            condense_water_cld,
            rh_limit);

    const Layout& layout = solver.get_layout();

    // Small grid with one ghost cell on each side, as MicroHH.
    const int itot = 8;
    const int jtot = 8;
    const int ktot = 8;
    const int gc = 1;

    const int icells = itot + 2*gc;
    const int jcells = jtot + 2*gc;
    const int kcells = ktot + 2*gc;

    Array_4d<TF> fields("fields", layout.n_fields, kcells, jcells, icells);
    Array_4d<TF> diameter_wet("diameter_wet", layout.n_bins, kcells, jcells, icells);
    Array_3d<TF> rv("rv", kcells, jcells, icells);
    Array_3d<TF> T("T", kcells, jcells, icells);
    Array_3d<TF> p("p", kcells, jcells, icells);
    Array_3d<TF> rs("rs", kcells, jcells, icells);
    Array_3d<TF> rho("rho", kcells, jcells, icells);

    const TF dt = 1._tf;
    const bool initialize = false;

    const bool do_coagulation = false;
    const bool update_coagulation_kernels = false;
    const bool do_condensation = false;
    const bool do_activation = false;
    const bool do_autoconversion = false;
    const bool do_update_distribution = false;
    const bool free_rh = true;

    const int istart = gc;
    const int iend = gc + itot;
    const int jstart = gc;
    const int jend = gc + jtot;
    const int kstart = gc;
    const int kend = gc + ktot;

    solver.exec(
            fields,
            diameter_wet,
            rv,
            T,
            p,
            rs,
            rho,
            dt,
            initialize,
            do_coagulation,
            update_coagulation_kernels,
            do_condensation,
            do_activation,
            do_autoconversion,
            do_update_distribution,
            free_rh,
            istart,
            iend,
            jstart,
            jend,
            kstart,
            kend);

    Kokkos::fence();

    std::printf("SALSA++ test: %s, %s, %d bins, %d species, %d fields, exec done\n",
            Default_exec::name(), sizeof(TF) == 4 ? "float" : "double",
            layout.n_bins, layout.n_species, layout.n_fields);

    return 0;
}
