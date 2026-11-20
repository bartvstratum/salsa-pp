/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * 0D driver: runs the solver at one grid point with the forcing and initial conditions of a reference case, with the
 * time stepping of the Fortran driver (fortran/driver/salsa_0d.f90), and writes output_cpp.nc in the same layout.
 * Reads the file NAMELIST (&salsa, &driver) in the current directory.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "salsa_solver.h"
#include "salsa_constants.h"
#include "salsa_sizedist.h"
#include "salsa_0d_io.h"
#include "salsa_0d_namelist.h"


namespace
{
    using namespace Salsa;
    using namespace Salsa_0d;

    // Saturation vapour pressure over liquid water (Pa): 6th-order polynomial in T - 273.15 K (Lowe, 1977), as the
    // Fortran driver.
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


    // Initial aerosol from lognormal modes, as the Fortran driver: number per bin for regimes 1a + 2a and for 2b,
    // mass = volume fraction * number * bin centre volume * density. Regime 1a holds only SO4 and OC, in the ratio of
    // their regime-a fractions.
    void set_initial_aerosol(
            const Array_4d_h& fields,
            const std::vector<TF>& number_a,
            const std::vector<TF>& dpg_a,
            const std::vector<TF>& sigmag_a,
            const std::vector<TF>& number_b,
            const std::vector<TF>& dpg_b,
            const std::vector<TF>& sigmag_b,
            const std::vector<TF>& vol_frac_a_in,
            const std::vector<TF>& vol_frac_b_in,
            const std::vector<std::string>& species_dry_name_in,
            const Layout& layout)
    {
        const int n_dry = layout.n_species - 1;

        // Volume fractions of the active dry species, by name.
        std::vector<TF> vol_frac_a(n_dry, 0._tf);
        std::vector<TF> vol_frac_b(n_dry, 0._tf);

        for (size_t j=0; j<species_dry_name_in.size(); ++j)
        {
            bool found = false;
            for (int s=0; s<n_dry; ++s)
            {
                if (layout.species[s] == species_dry_name_in[j])
                {
                    vol_frac_a[s] = vol_frac_a_in[j];
                    vol_frac_b[s] = vol_frac_b_in[j];
                    found = true;
                }
            }

            if (!found && (vol_frac_a_in[j] != 0._tf || vol_frac_b_in[j] != 0._tf))
                throw std::runtime_error("Initial volume fraction given for " + species_dry_name_in[j] + ", which is not in listspec");
        }

        const std::vector<TF> number_bin_a = lognormal_distribution(
                number_a, dpg_a, sigmag_a, layout, layout.aer_1a_start, layout.aer_2a_end);
        const std::vector<TF> number_bin_b = lognormal_distribution(
                number_b, dpg_b, sigmag_b, layout, layout.aer_2b_start, layout.aer_2b_end);

        int iso4 = -1;
        int ioc = -1;
        for (int s=0; s<n_dry; ++s)
        {
            if (layout.species[s] == "SO4")
                iso4 = s;
            if (layout.species[s] == "OC")
                ioc = s;
        }

        TF frac_oc_1a;
        if (iso4 >= 0 && ioc >= 0)
            frac_oc_1a = vol_frac_a[ioc]/(vol_frac_a[ioc] + vol_frac_a[iso4]);
        else if (ioc >= 0)
            frac_oc_1a = 1._tf;
        else if (iso4 >= 0)
            frac_oc_1a = 0._tf;
        else
            throw std::runtime_error("Either OC or SO4 must be active for aerosol regime 1a");

        for (int b=layout.aer_1a_start; b<layout.aer_2b_end; ++b)
        {
            const TF number = number_bin_a[b] + number_bin_b[b];
            const TF d = layout.diameter_mid[b];
            const TF core = pi6*(d*d*d);

            fields(layout.number_field(b), 0, 0, 0) = number;

            if (b < layout.aer_1a_end)
            {
                if (iso4 >= 0)
                    fields(layout.mass_field(b, iso4), 0, 0, 0) =
                            std::max(0._tf, 1._tf - frac_oc_1a)*number*core*layout.density[iso4];
                if (ioc >= 0)
                    fields(layout.mass_field(b, ioc), 0, 0, 0) =
                            std::max(0._tf, frac_oc_1a)*number*core*layout.density[ioc];
            }
            else
            {
                const std::vector<TF>& vol_frac = (b < layout.aer_2a_end) ? vol_frac_a : vol_frac_b;
                for (int s=0; s<n_dry; ++s)
                    fields(layout.mass_field(b, s), 0, 0, 0) =
                            std::max(0._tf, vol_frac[s])*number*core*layout.density[s];
            }
        }
    }


    std::string read_text_file(const std::string& file_name)
    {
        std::ifstream file(file_name);
        std::stringstream text;
        text << file.rdbuf();
        return text.str();
    }


    void run()
    {
        const std::string namelist_file = "NAMELIST";
        const std::string output_file = "output_cpp.nc";

        const Namelist_group salsa = read_namelist_group(namelist_file, "salsa");
        const Namelist_group driver = read_namelist_group(namelist_file, "driver");

        // Items of &salsa that this driver supports; anything else would silently differ from the Fortran.
        const std::set<std::string> supported = {
                "reglim", "nbin", "bloprc%nbins", "bloprc%dlo", "bloprc%vol_ratio", "nspec_dry", "listspec",
                "lscoag%switch", "lscoag%mode", "lscoag%delay", "lscnd%switch", "lscnd%delay",
                "lsauto%switch", "lsauto%mode", "lsauto%delay", "lsactiv%switch", "lsactiv%mode", "lsactiv%delay",
                "lscgaa", "lscgcc", "lscgca", "lscgpc", "lscgpa", "lscgpp", "cgintvl",
                "lscndgas", "lscndh2oae", "lscndh2ocl", "lsdistupdate", "lsfreerh%switch", "lsfreerh%delay", "rhlim"};

        for (const auto& item : salsa)
            if (supported.find(item.first) == supported.end())
                throw std::runtime_error("Namelist item " + item.first + " in &salsa is not supported by the C++ driver");

        // &driver.
        const std::string case_name = get_item<std::string>(driver, "case_name");
        const std::string input_file = get_item<std::string>(driver, "input_file", "input.nc");
        const TF dt = get_item<TF>(driver, "dt", 1._tf);
        const int nsteps = get_item<int>(driver, "nsteps", 400);

        // &salsa, with the defaults of mo_submctl.f90 where the reference cases do not set an item.
        const std::vector<TF> reglim = get_list<TF>(salsa, "reglim");
        const std::vector<int> nbin = get_list<int>(salsa, "nbin");
        if (reglim.size() != 4 || nbin.size() != 2)
            throw std::runtime_error("reglim needs 4 values and nbin 2");

        const std::array<TF, 4> regime_limits = {reglim[0], reglim[1], reglim[2], reglim[3]};
        const std::array<int, 2> n_bins_regime = {nbin[0], nbin[1]};
        const int n_bins_rain = get_item<int>(salsa, "bloprc%nbins");
        const TF rain_diameter_low = get_item<TF>(salsa, "bloprc%dlo");
        const TF rain_volume_ratio = get_item<TF>(salsa, "bloprc%vol_ratio");

        std::vector<std::string> dry_species;
        for (const std::string& name : get_list<std::string>(salsa, "listspec"))
            if (!name.empty())
                dry_species.push_back(name);

        if (static_cast<int>(dry_species.size()) != get_item<int>(salsa, "nspec_dry"))
            throw std::runtime_error("nspec_dry does not match the number of names in listspec");

        const bool coagulation_on = get_item<bool>(salsa, "lscoag%switch", false);
        const int coagulation_mode = get_item<int>(salsa, "lscoag%mode", 1);
        const TF coagulation_delay = get_item<TF>(salsa, "lscoag%delay", 0._tf);
        const TF coagulation_interval = get_item<TF>(salsa, "cgintvl", 10._tf);

        const bool condensation_on = get_item<bool>(salsa, "lscnd%switch", false);
        const TF condensation_delay = get_item<TF>(salsa, "lscnd%delay", 0._tf);

        const bool autoconversion_on = get_item<bool>(salsa, "lsauto%switch", false);
        const int autoconversion_mode = get_item<int>(salsa, "lsauto%mode", 1);
        const TF autoconversion_delay = get_item<TF>(salsa, "lsauto%delay", 0._tf);

        const bool activation_on = get_item<bool>(salsa, "lsactiv%switch", false);
        const int activation_mode = get_item<int>(salsa, "lsactiv%mode", 1);
        const TF activation_delay = get_item<TF>(salsa, "lsactiv%delay", 0._tf);

        const bool coagulate_aer = get_item<bool>(salsa, "lscgaa", true);
        const bool coagulate_cld = get_item<bool>(salsa, "lscgcc", true);
        const bool cld_collects_aer = get_item<bool>(salsa, "lscgca", true);
        const bool rain_collects_cld = get_item<bool>(salsa, "lscgpc", true);
        const bool rain_collects_aer = get_item<bool>(salsa, "lscgpa", true);
        const bool coagulate_rain = get_item<bool>(salsa, "lscgpp", true);

        const bool condense_gases = get_item<bool>(salsa, "lscndgas", false);
        const bool condense_water_aer = get_item<bool>(salsa, "lscndh2oae", true);
        const bool condense_water_cld = get_item<bool>(salsa, "lscndh2ocl", true);

        const bool do_update_distribution = get_item<bool>(salsa, "lsdistupdate", true);
        const bool free_rh_on = get_item<bool>(salsa, "lsfreerh%switch", false);
        const TF free_rh_delay = get_item<TF>(salsa, "lsfreerh%delay", 0._tf);
        const TF rh_limit = get_item<TF>(salsa, "rhlim", 1.20_tf);

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

        // Forcing and initial conditions.
        const std::vector<double> time_in = read_variable(input_file, "time");
        const std::vector<double> T_in = read_variable(input_file, "temperature");
        const std::vector<double> p_in = read_variable(input_file, "pressure");
        const std::vector<double> rh_in = read_variable(input_file, "relative_humidity");
        const std::vector<double> gas_in = read_variable(input_file, "gas");
        const std::vector<std::string> gas_name_in = read_names(input_file, "gas_name");

        auto read_tf = [&](const std::string& name)
        {
            const std::vector<double> values = read_variable(input_file, name);
            return std::vector<TF>(values.begin(), values.end());
        };

        const std::vector<TF> number_a = read_tf("number_a");
        const std::vector<TF> dpg_a = read_tf("dpg_a");
        const std::vector<TF> sigmag_a = read_tf("sigmag_a");
        const std::vector<TF> number_b = read_tf("number_b");
        const std::vector<TF> dpg_b = read_tf("dpg_b");
        const std::vector<TF> sigmag_b = read_tf("sigmag_b");
        const std::vector<TF> vol_frac_a_in = read_tf("vol_frac_a");
        const std::vector<TF> vol_frac_b_in = read_tf("vol_frac_b");
        const std::vector<std::string> species_dry_name_in = read_names(input_file, "species_dry_name");

        const std::vector<std::string> gas_names = {"H2SO4", "HNO3", "NH3", "OCNV", "OCSV"};

        if (static_cast<int>(time_in.size()) < nsteps + 1)
            throw std::runtime_error("The input file has fewer than nsteps + 1 forcing times");
        if (gas_name_in != gas_names)
            throw std::runtime_error("Unexpected gases or gas order in the input file");

        // One grid point.
        const int istart = 0;
        const int iend = 1;
        const int jstart = 0;
        const int jend = 1;
        const int kstart = 0;
        const int kend = 1;

        Array_4d<TF> fields("fields", layout.n_fields, kend, jend, iend);
        Array_4d<TF> diameter_wet("diameter_wet", layout.n_bins, kend, jend, iend);
        Array_3d<TF> rv("rv", kend, jend, iend);
        Array_3d<TF> T("T", kend, jend, iend);
        Array_3d<TF> p("p", kend, jend, iend);
        Array_3d<TF> rs("rs", kend, jend, iend);
        Array_3d<TF> rho("rho", kend, jend, iend);

        auto fields_h = Kokkos::create_mirror_view(fields);
        auto rv_h = Kokkos::create_mirror_view(rv);
        auto T_h = Kokkos::create_mirror_view(T);
        auto p_h = Kokkos::create_mirror_view(p);
        auto rs_h = Kokkos::create_mirror_view(rs);
        auto rho_h = Kokkos::create_mirror_view(rho);

        // Temperature and pressure from the forcing; air density and saturation mixing ratio as the Fortran driver.
        auto set_thermodynamics = [&](const int row)
        {
            T_h(0, 0, 0) = static_cast<TF>(T_in[row]);
            p_h(0, 0, 0) = static_cast<TF>(p_in[row]);
            rho_h(0, 0, 0) = mair*p_h(0, 0, 0)/(rg*T_h(0, 0, 0));
            rs_h(0, 0, 0) = mwa*saturation_vapour_pressure(T_h(0, 0, 0))/(mair*p_h(0, 0, 0));

            Kokkos::deep_copy(T, T_h);
            Kokkos::deep_copy(p, p_h);
            Kokkos::deep_copy(rho, rho_h);
            Kokkos::deep_copy(rs, rs_h);
        };

        // Which processes run at this time, as set_SALSA_runtime in mo_salsa_driver.f90.
        auto exec = [&](const TF time, const int istp)
        {
            const bool do_coagulation = coagulation_on && time >= coagulation_delay;
            const bool update_coagulation_kernels = do_coagulation && coagulation_mode == 2 &&
                    (std::fmod(time, coagulation_interval) < dt || istp <= 1);
            const bool do_condensation = condensation_on && time >= condensation_delay;
            const bool do_activation = activation_on && time >= activation_delay;
            const bool do_autoconversion = autoconversion_on && time >= autoconversion_delay;
            const bool free_rh = free_rh_on && time > free_rh_delay;

            solver.exec(
                    fields,
                    diameter_wet,
                    rv,
                    T,
                    p,
                    rs,
                    rho,
                    dt,
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

        Output output(output_file, case_name, input_file, read_text_file(namelist_file), layout);

        auto write_output = [&](const TF time)
        {
            Kokkos::deep_copy(fields_h, fields);
            Kokkos::deep_copy(rv_h, rv);

            output.write(
                    fields_h,
                    time,
                    T_h(0, 0, 0),
                    p_h(0, 0, 0),
                    rho_h(0, 0, 0),
                    rv_h(0, 0, 0),
                    rs_h(0, 0, 0));
        };

        // Initial state.
        set_thermodynamics(0);
        rv_h(0, 0, 0) = static_cast<TF>(rh_in[0])*rs_h(0, 0, 0);
        for (int g=0; g<Layout::n_gases; ++g)
            fields_h(layout.gas_field(g), 0, 0, 0) = static_cast<TF>(gas_in[g]);

        set_initial_aerosol(
                fields_h,
                number_a,
                dpg_a,
                sigmag_a,
                number_b,
                dpg_b,
                sigmag_b,
                vol_frac_a_in,
                vol_frac_b_in,
                species_dry_name_in,
                layout);

        // The bin centre diameter, as UCLALES-SALSA initialises the wet diameter.
        auto diameter_wet_h = Kokkos::create_mirror_view(diameter_wet);
        for (int b=0; b<layout.n_bins; ++b)
            diameter_wet_h(b, 0, 0, 0) = layout.diameter_mid[b];

        Kokkos::deep_copy(rv, rv_h);
        Kokkos::deep_copy(fields, fields_h);
        Kokkos::deep_copy(diameter_wet, diameter_wet_h);

        // Initialisation call at time 0, as the Fortran driver.
        const TF time_start = 0._tf;
        const int istp_start = 0;
        exec(time_start, istp_start);
        write_output(time_start);

        // Step istp sees the forcing at its end time istp*dt, and gets the time at its start.
        for (int istp=1; istp<=nsteps; ++istp)
        {
            if (std::abs(time_in[istp] - istp*dt) > 1.e-6*dt)
                throw std::runtime_error("Input time does not match istp*dt");

            set_thermodynamics(istp);
            exec((istp-1)*dt, istp);
            write_output(istp*dt);
        }

        std::printf("salsa_0d: case %s, %d steps, output in %s\n", case_name.c_str(), nsteps, output_file.c_str());
    }
}


int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);

    try
    {
        run();
    }
    catch (const std::exception& e)
    {
        std::printf("salsa_0d: %s\n", e.what());
        return 1;
    }

    return 0;
}
