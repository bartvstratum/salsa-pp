/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Benchmark: times the solver per process on a 3D grid, with the bins and aerosol of case C without rain.
 *
 * Usage: salsa_bench itot jtot ktot dt [n_calls]
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "salsa_solver.h"
#include "salsa_driver_time_integration.h"
#include "salsa_constants.h"
#include "salsa_sizedist.h"


namespace
{
    using namespace Salsa;

    // Saturation vapour pressure (Pa), as the 0D driver.
    TF saturation_vapour_pressure(const TF T)
    {
        const TF a0 = 6.107799961_tf;
        const TF a1 = 4.436518521e-1_tf;
        const TF a2 = 1.428945805e-2_tf;
        const TF a3 = 2.650648471e-4_tf;
        const TF a4 = 3.031240396e-6_tf;
        const TF a5 = 2.034080948e-8_tf;
        const TF a6 = 6.136820929e-11_tf;

        const TF tc = T - 273.15_tf;
        const TF tc2 = tc*tc;
        const TF tc3 = tc2*tc;
        const TF tc4 = tc2*tc2;
        const TF tc5 = tc4*tc;
        const TF tc6 = tc3*tc3;

        return 100._tf*(a0 + a1*tc + a2*tc2 + a3*tc3 + a4*tc4 + a5*tc5 + a6*tc6);
    }


    // Aerosol of case C at one grid point (per kg of air).
    std::vector<TF> initial_aerosol(const Layout& layout)
    {
        const std::vector<TF> number_a = {1.e8_tf, 1.e6_tf};
        const std::vector<TF> dpg_a = {0.15e-6_tf, 1.0e-6_tf};
        const std::vector<TF> sigmag_a = {2.0_tf, 1.5_tf};
        const std::vector<TF> number_b = {1.e7_tf};
        const std::vector<TF> dpg_b = {0.5e-6_tf};
        const std::vector<TF> sigmag_b = {1.7_tf};

        // In the order of the dry species of the layout: SO4, OC, BC, DU, SS, NO, NH.
        const std::vector<TF> vol_frac_a = {0.30_tf, 0.25_tf, 0.05_tf, 0.05_tf, 0.15_tf, 0.10_tf, 0.10_tf};
        const std::vector<TF> vol_frac_b = {0.05_tf, 0.05_tf, 0.40_tf, 0.40_tf, 0.04_tf, 0.03_tf, 0.03_tf};

        const int iso4 = 0;
        const int ioc = 1;
        const TF frac_oc_1a = vol_frac_a[ioc]/(vol_frac_a[ioc] + vol_frac_a[iso4]);

        const std::vector<TF> number_bin_a = lognormal_distribution(
                number_a, dpg_a, sigmag_a, layout, layout.aer_1a_start, layout.aer_2a_end);
        const std::vector<TF> number_bin_b = lognormal_distribution(
                number_b, dpg_b, sigmag_b, layout, layout.aer_2b_start, layout.aer_2b_end);

        const int n_dry = layout.n_species - 1;
        std::vector<TF> values(layout.n_fields, 0._tf);

        for (int b=layout.aer_1a_start; b<layout.aer_2b_end; ++b)
        {
            const TF number = number_bin_a[b] + number_bin_b[b];
            const TF d = layout.diameter_mid[b];
            const TF core = pi6*(d*d*d);

            values[layout.number_field(b)] = number;

            if (b < layout.aer_1a_end)
            {
                values[layout.mass_field(b, iso4)] = (1._tf - frac_oc_1a)*number*core*layout.density[iso4];
                values[layout.mass_field(b, ioc)] = frac_oc_1a*number*core*layout.density[ioc];
            }
            else
            {
                const std::vector<TF>& vol_frac = (b < layout.aer_2a_end) ? vol_frac_a : vol_frac_b;
                for (int s=0; s<n_dry; ++s)
                    values[layout.mass_field(b, s)] = vol_frac[s]*number*core*layout.density[s];
            }
        }

        // HNO3 and NH3 as case C (molecules per kg).
        const TF gas_hno3 = 8.4682e-15_tf*avog*1.e6_tf;
        values[layout.gas_field(Gas::hno3)] = gas_hno3;
        values[layout.gas_field(Gas::nh3)] = gas_hno3;

        return values;
    }


    void run(int argc, char* argv[])
    {
        if (argc != 5 && argc != 6)
            throw std::runtime_error("usage: salsa_bench itot jtot ktot dt [n_calls]");

        const int itot = std::stoi(argv[1]);
        const int jtot = std::stoi(argv[2]);
        const int ktot = std::stoi(argv[3]);
        const TF dt = static_cast<TF>(std::stod(argv[4]));
        const int n_calls = (argc == 6) ? std::stoi(argv[5]) : 10;

        // Bins and species of case C, without rain (LES target).
        const std::array<TF, 4> regime_limits = {3.e-9_tf, 5.e-8_tf, 7.e-7_tf, 1.e-5_tf};
        const std::array<int, 2> n_bins_regime = {3, 7};
        const int n_bins_rain = 0;
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

        // Grid with one ghost cell on each side, as MicroHH.
        const int gc = 1;

        const int icells = itot + 2*gc;
        const int jcells = jtot + 2*gc;
        const int kcells = ktot + 2*gc;

        const int istart = gc;
        const int iend = gc + itot;
        const int jstart = gc;
        const int jend = gc + jtot;
        const int kstart = gc;
        const int kend = gc + ktot;

        Array_4d<TF> fields("fields", layout.n_fields, kcells, jcells, icells);
        Array_4d<TF> fields_tend("fields_tend", layout.n_fields, kcells, jcells, icells);
        Array_4d<TF> diameter_wet("diameter_wet", layout.n_bins, kcells, jcells, icells);
        Array_3d<TF> rv("rv", kcells, jcells, icells);
        Array_3d<TF> rv_tend("rv_tend", kcells, jcells, icells);
        Array_3d<TF> T("T", kcells, jcells, icells);
        Array_3d<TF> p("p", kcells, jcells, icells);
        Array_3d<TF> rs("rs", kcells, jcells, icells);
        Array_3d<TF> rho("rho", kcells, jcells, icells);

        auto fields_h = Kokkos::create_mirror_view(fields);
        auto diameter_wet_h = Kokkos::create_mirror_view(diameter_wet);
        auto rv_h = Kokkos::create_mirror_view(rv);
        auto T_h = Kokkos::create_mirror_view(T);
        auto p_h = Kokkos::create_mirror_view(p);
        auto rs_h = Kokkos::create_mirror_view(rs);
        auto rho_h = Kokkos::create_mirror_view(rho);

        // 3 km domain: RH 95 % below 1 km, 100.2 % up to 2 km, 60 % above, plus noise.
        const TF domain_height = 3000._tf;
        const TF T_surface = 290._tf;
        const TF lapse_rate = 0.0065_tf;
        const TF p_surface = 1.e5_tf;
        const TF scale_height = 8000._tf;
        const TF cloud_base = 1000._tf;
        const TF cloud_top = 2000._tf;
        const TF rh_below = 0.95_tf;
        const TF rh_cloud = 1.002_tf;
        const TF rh_above = 0.60_tf;
        const TF rh_perturbation = 0.005_tf;

        std::mt19937 random_engine(1);
        std::uniform_real_distribution<double> random_uniform(-1., 1.);

        const std::vector<TF> aerosol = initial_aerosol(layout);

        for (int k=kstart; k<kend; ++k)
        {
            const TF z = (k - kstart + 0.5_tf)/ktot*domain_height;
            const TF rh_z = (z < cloud_base) ? rh_below : (z < cloud_top) ? rh_cloud : rh_above;

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    T_h(k, j, i) = T_surface - lapse_rate*z;
                    p_h(k, j, i) = p_surface*std::exp(-z/scale_height);
                    rho_h(k, j, i) = mair*p_h(k, j, i)/(rg*T_h(k, j, i));
                    rs_h(k, j, i) = mwa*saturation_vapour_pressure(T_h(k, j, i))/(mair*p_h(k, j, i));

                    const TF rh = rh_z + rh_perturbation*static_cast<TF>(random_uniform(random_engine));
                    rv_h(k, j, i) = rh*rs_h(k, j, i);

                    for (int n=0; n<layout.n_fields; ++n)
                        fields_h(n, k, j, i) = aerosol[n];

                    // The bin centre diameter, as UCLALES-SALSA initialises the wet diameter.
                    for (int b=0; b<layout.n_bins; ++b)
                        diameter_wet_h(b, k, j, i) = layout.diameter_mid[b];
                }
        }

        Kokkos::deep_copy(fields, fields_h);
        Kokkos::deep_copy(diameter_wet, diameter_wet_h);
        Kokkos::deep_copy(rv, rv_h);
        Kokkos::deep_copy(T, T_h);
        Kokkos::deep_copy(p, p_h);
        Kokkos::deep_copy(rs, rs_h);
        Kokkos::deep_copy(rho, rho_h);

        const bool update_coagulation_kernels = false;
        const bool do_coagulation = false;
        const bool do_autoconversion = false;
        const bool free_rh = true;

        auto exec = [&](
                const bool initialize,
                const bool do_condensation,
                const bool do_activation,
                const bool do_update_distribution)
        {
            solver.exec(
                    fields_tend,
                    rv_tend,
                    diameter_wet,
                    fields,
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
        };

        // The fields evolve during the spin-up only, as in the 0D driver.
        auto integrate = [&](const bool initialize)
        {
            if (!initialize)
                constrain_tendencies(fields, fields_tend, dt, istart, iend, jstart, jend, kstart, kend);

            integrate_forward_euler(fields, fields_tend, rv, rv_tend, dt, initialize, istart, iend, jstart, jend, kstart, kend);
        };

        // Spin-up: initialisation, then 10 s of all processes to fill the cloud bins.
        const TF spin_up_time = 10._tf;
        const int n_spin_up = std::max(1, static_cast<int>(std::lround(spin_up_time/dt)));

        const bool initialize = true;
        const bool all_off = false;
        const bool all_on = true;

        exec(initialize, all_off, all_off, all_off);
        integrate(initialize);
        for (int n=0; n<n_spin_up; ++n)
        {
            exec(!initialize, all_on, all_on, all_on);
            integrate(!initialize);
        }

        // Start state of every timed call.
        Array_4d<TF> fields_start("fields_start", layout.n_fields, kcells, jcells, icells);
        Array_4d<TF> diameter_wet_start("diameter_wet_start", layout.n_bins, kcells, jcells, icells);
        Array_3d<TF> rv_start("rv_start", kcells, jcells, icells);

        Kokkos::deep_copy(fields_start, fields);
        Kokkos::deep_copy(diameter_wet_start, diameter_wet);
        Kokkos::deep_copy(rv_start, rv);

        const std::vector<std::string> process_names = {"none", "condensation", "activation", "distr_update", "all"};
        const std::vector<bool> condensation_on = {false, true, false, false, true};
        const std::vector<bool> activation_on = {false, false, true, false, true};
        const std::vector<bool> update_distribution_on = {false, false, false, true, true};

        const double n_points = static_cast<double>(itot)*jtot*ktot;

        std::printf("SALSA++ benchmark: %s, %s, %d x %d x %d = %.0f points, dt = %g s, %d fields, %d calls\n",
                Default_exec::name(), sizeof(TF) == 4 ? "float" : "double",
                itot, jtot, ktot, n_points, static_cast<double>(dt), layout.n_fields, n_calls);
        std::printf("%-14s %12s %16s %18s\n", "process", "ms/call", "ns/point/call", "ms/step (3 calls)");

        for (size_t n=0; n<process_names.size(); ++n)
        {
            // One warm-up call, then n_calls timed calls, each from the same start state.
            double time_total = 0.;

            for (int c=0; c<n_calls+1; ++c)
            {
                Kokkos::deep_copy(fields, fields_start);
                Kokkos::deep_copy(diameter_wet, diameter_wet_start);
                Kokkos::deep_copy(rv, rv_start);
                Kokkos::deep_copy(fields_tend, 0._tf);
                Kokkos::deep_copy(rv_tend, 0._tf);
                Kokkos::fence();

                Kokkos::Timer timer;
                exec(!initialize, condensation_on[n], activation_on[n], update_distribution_on[n]);
                Kokkos::fence();

                if (c > 0)
                    time_total += timer.seconds();
            }

            const double time_call = time_total/n_calls;

            std::printf("%-14s %12.3f %16.1f %18.3f\n",
                    process_names[n].c_str(), 1.e3*time_call, 1.e9*time_call/n_points, 3.e3*time_call);
        }
    }
}


int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);

    try
    {
        run(argc, argv);
    }
    catch (const std::exception& e)
    {
        std::printf("salsa_bench: %s\n", e.what());
        return 1;
    }

    return 0;
}
