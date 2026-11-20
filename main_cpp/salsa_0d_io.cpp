/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * NetCDF input and output of the 0D driver. The output has the layout of the Fortran driver's output.
 */

#include <cmath>
#include <ctime>
#include <stdexcept>

#include <netcdf.h>

#include "salsa_0d_io.h"
#include "salsa_constants.h"
#include "salsa_version.h"


namespace
{
    constexpr int name_strlen = 8;

    const std::array<std::string, 3> categories = {"aerosol", "cloud", "rain"};
    const std::array<std::string, Salsa::Gas::count> gas_names = {"H2SO4", "HNO3", "NH3", "OCNV", "OCSV"};


    void check(const int status, const std::string& what)
    {
        if (status != NC_NOERR)
            throw std::runtime_error("NetCDF error (" + what + "): " + nc_strerror(status));
    }


    // Without the trailing spaces and NULs of a Fortran or NumPy character array.
    std::string strip(const std::string& name)
    {
        const size_t last = name.find_last_not_of(std::string(" \0", 2));
        return last == std::string::npos ? "" : name.substr(0, last + 1);
    }


    void put_text(const int ncid, const int varid, const std::string& name, const std::string& value)
    {
        check(nc_put_att_text(ncid, varid, name.c_str(), value.size(), value.c_str()), name);
    }


    int def_var(
            const int ncid,
            const std::string& name,
            const std::vector<int>& dimids,
            const std::string& units,
            const std::string& long_name,
            const std::string& standard_name = "")
    {
        int varid;
        check(nc_def_var(ncid, name.c_str(), NC_DOUBLE, static_cast<int>(dimids.size()), dimids.data(), &varid), name);

        put_text(ncid, varid, "units", units);
        put_text(ncid, varid, "long_name", long_name);
        if (!standard_name.empty())
            put_text(ncid, varid, "standard_name", standard_name);

        return varid;
    }


    void put_name(const int ncid, const int varid, const int index, const std::string& name)
    {
        std::string padded = name;
        padded.resize(name_strlen, ' ');

        const std::array<size_t, 2> start = {static_cast<size_t>(index), 0};
        const std::array<size_t, 2> count = {1, name_strlen};
        check(nc_put_vara_text(ncid, varid, start.data(), count.data(), padded.data()), name);
    }
}


namespace Salsa_0d
{
    std::vector<double> read_variable(const std::string& file_name, const std::string& name)
    {
        int ncid;
        int varid;
        int ndims;
        int dimid;
        size_t n;

        check(nc_open(file_name.c_str(), NC_NOWRITE, &ncid), file_name);
        check(nc_inq_varid(ncid, name.c_str(), &varid), name);
        check(nc_inq_varndims(ncid, varid, &ndims), name);
        if (ndims != 1)
            throw std::runtime_error("Input variable " + name + " is not 1D");

        check(nc_inq_vardimid(ncid, varid, &dimid), name);
        check(nc_inq_dimlen(ncid, dimid, &n), name);

        std::vector<double> values(n);
        check(nc_get_var_double(ncid, varid, values.data()), name);
        check(nc_close(ncid), file_name);

        return values;
    }


    std::vector<std::string> read_names(const std::string& file_name, const std::string& name)
    {
        int ncid;
        int varid;
        int ndims;
        std::array<int, 2> dimids;
        size_t n;
        size_t length;

        check(nc_open(file_name.c_str(), NC_NOWRITE, &ncid), file_name);
        check(nc_inq_varid(ncid, name.c_str(), &varid), name);
        check(nc_inq_varndims(ncid, varid, &ndims), name);
        if (ndims != 2)
            throw std::runtime_error("Input variable " + name + " is not a list of names");

        check(nc_inq_vardimid(ncid, varid, dimids.data()), name);
        check(nc_inq_dimlen(ncid, dimids[0], &n), name);
        check(nc_inq_dimlen(ncid, dimids[1], &length), name);

        std::vector<std::string> names;
        for (size_t i=0; i<n; ++i)
        {
            std::string text(length, ' ');
            const std::array<size_t, 2> start = {i, 0};
            const std::array<size_t, 2> count = {1, length};
            check(nc_get_vara_text(ncid, varid, start.data(), count.data(), text.data()), name);
            names.push_back(strip(text));
        }

        check(nc_close(ncid), file_name);

        return names;
    }


    Output::Output(
            const std::string& file_name,
            const std::string& case_name,
            const std::string& input_file,
            const std::string& namelist_text,
            const Salsa::Layout& layout) :
        layout(layout),
        category_start{layout.aer_1a_start, layout.cld_a_start, layout.rain_start},
        category_end{layout.aer_2b_end, layout.cld_b_end, layout.rain_end},
        record(0)
    {
        check(nc_create(file_name.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &ncid), file_name);

        char creation_time[20];
        const std::time_t now = std::time(nullptr);
        std::strftime(creation_time, sizeof(creation_time), "%Y-%m-%dT%H:%M:%S", std::localtime(&now));

        put_text(ncid, NC_GLOBAL, "title", "SALSA 0D reference test, case " + case_name);
        put_text(ncid, NC_GLOBAL, "case_name", case_name);
        put_text(ncid, NC_GLOBAL, "code", "SALSA++ (C++) 0D driver");
        put_text(ncid, NC_GLOBAL, "salsapp_git_commit", salsapp_git_commit);
        put_text(ncid, NC_GLOBAL, "upstream_salsa_commit",
                "https://github.com/UCLALES-SALSA/UCLALES-SALSA 8c779581118f97b6092efe6da43680ffa4e1cefa");
        put_text(ncid, NC_GLOBAL, "creation_time", creation_time);
        put_text(ncid, NC_GLOBAL, "input_file", input_file);
        put_text(ncid, NC_GLOBAL, "namelist", namelist_text);
        put_text(ncid, NC_GLOBAL, "Conventions", "CF-1.8");

        const int n_bounds = 2;

        int dimid_time;
        int dimid_species;
        int dimid_gas;
        int dimid_nv;
        int dimid_strlen;
        std::array<int, 3> dimid_bin;

        check(nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time), "time");
        for (int c=0; c<3; ++c)
        {
            const std::string name = "bin_" + categories[c];
            check(nc_def_dim(ncid, name.c_str(), category_end[c] - category_start[c], &dimid_bin[c]), name);
        }
        check(nc_def_dim(ncid, "species", layout.n_species, &dimid_species), "species");
        check(nc_def_dim(ncid, "gas", Salsa::Gas::count, &dimid_gas), "gas");
        check(nc_def_dim(ncid, "nv", n_bounds, &dimid_nv), "nv");
        check(nc_def_dim(ncid, "name_strlen", name_strlen, &dimid_strlen), "name_strlen");

        varid_time = def_var(ncid, "time", {dimid_time}, "s", "time since start");

        // Bin coordinates: dry diameter for aerosol and cloud (cloud bins are parallel to aerosol bins), drop diameter for rain.
        std::array<int, 3> varid_bin;
        std::array<int, 3> varid_bounds;

        for (int c=0; c<3; ++c)
        {
            const std::string name = "bin_" + categories[c];
            const std::string long_name = (c < 2) ?
                    categories[c] + " bin centre dry diameter (SALSA dmid)" :
                    "rain bin centre drop diameter (SALSA dmid)";

            varid_bin[c] = def_var(ncid, name, {dimid_bin[c]}, "m", long_name);
            put_text(ncid, varid_bin[c], "bounds", name + "_bounds");
            varid_bounds[c] = def_var(ncid, name + "_bounds", {dimid_bin[c], dimid_nv}, "m", categories[c] + " bin diameter limits");
        }

        int varid_species;
        int varid_gas_name;
        const std::array<int, 2> dimids_species = {dimid_species, dimid_strlen};
        const std::array<int, 2> dimids_gas_name = {dimid_gas, dimid_strlen};

        check(nc_def_var(ncid, "species", NC_CHAR, 2, dimids_species.data(), &varid_species), "species");
        put_text(ncid, varid_species, "long_name", "compound names, in the order of the species dimension");
        check(nc_def_var(ncid, "gas_name", NC_CHAR, 2, dimids_gas_name.data(), &varid_gas_name), "gas_name");
        put_text(ncid, varid_gas_name, "long_name", "gas names, in the order of the gas dimension");

        varid_T = def_var(ncid, "temperature", {dimid_time}, "K", "air temperature", "air_temperature");
        varid_p = def_var(ncid, "pressure", {dimid_time}, "Pa", "air pressure", "air_pressure");
        varid_rho = def_var(ncid, "air_density", {dimid_time}, "kg m-3", "dry air density, M_air p / (R T)");
        varid_rv = def_var(ncid, "qv", {dimid_time}, "kg kg-1", "water vapour mixing ratio", "humidity_mixing_ratio");
        varid_rs = def_var(ncid, "qsat", {dimid_time}, "kg kg-1",
                "saturation water vapour mixing ratio, M_w e_s(T) / (M_air p), e_s from a 6th-order polynomial");
        varid_supersaturation = def_var(ncid, "supersaturation", {dimid_time}, "1",
                "supersaturation with respect to liquid water, qv/qsat - 1");
        varid_gas = def_var(ncid, "gas", {dimid_time, dimid_gas}, "kg-1", "gas concentration, molecules per kg of air");

        const double fill = NC_FILL_DOUBLE;

        for (int c=0; c<3; ++c)
        {
            const std::string& category = categories[c];

            varid_number[c] = def_var(ncid, "number_" + category, {dimid_time, dimid_bin[c]}, "kg-1",
                    category + " number per kg of air");
            varid_mass[c] = def_var(ncid, "mass_" + category, {dimid_time, dimid_species, dimid_bin[c]}, "kg kg-1",
                    category + " mass per compound per kg of air");
            varid_diameter_wet[c] = def_var(ncid, "diameter_wet_" + category, {dimid_time, dimid_bin[c]}, "m",
                    category + " wet diameter: volume-equivalent diameter of the mean particle, all compounds");
            varid_diameter_dry[c] = def_var(ncid, "diameter_dry_" + category, {dimid_time, dimid_bin[c]}, "m",
                    category + " dry diameter: volume-equivalent diameter of the mean particle, dry compounds");

            check(nc_put_att_double(ncid, varid_diameter_wet[c], "_FillValue", NC_DOUBLE, 1, &fill), "_FillValue");
            check(nc_put_att_double(ncid, varid_diameter_dry[c], "_FillValue", NC_DOUBLE, 1, &fill), "_FillValue");
            put_text(ncid, varid_diameter_wet[c], "description",
                    "Fill value where the bin is empty (number below the SALSA limit nlim, or prlim for rain).");
        }

        check(nc_enddef(ncid), "enddef");

        for (int c=0; c<3; ++c)
        {
            const size_t n_bins = category_end[c] - category_start[c];

            Kokkos::View<double*, Kokkos::HostSpace> diameter_mid("diameter_mid", n_bins);
            Kokkos::View<double**, Kokkos::HostSpace> bounds("bounds", n_bins, n_bounds);

            for (size_t b=0; b<n_bins; ++b)
            {
                const int bin = category_start[c] + b;
                diameter_mid(b) = layout.diameter_mid[bin];
                bounds(b, 0) = std::pow(layout.volume_low[bin]/Salsa::pi6, 1./3.);
                bounds(b, 1) = std::pow(layout.volume_high[bin]/Salsa::pi6, 1./3.);
            }

            check(nc_put_var_double(ncid, varid_bin[c], diameter_mid.data()), "bin");
            check(nc_put_var_double(ncid, varid_bounds[c], bounds.data()), "bounds");
        }

        for (int s=0; s<layout.n_species; ++s)
            put_name(ncid, varid_species, s, layout.species[s]);

        for (int g=0; g<Salsa::Gas::count; ++g)
            put_name(ncid, varid_gas_name, g, gas_names[g]);
    }


    Output::~Output()
    {
        nc_close(ncid);
    }


    void Output::write(
            const Array_4d_h& fields,
            const double time,
            const double T,
            const double p,
            const double rho,
            const double rv,
            const double rs)
    {
        const size_t rec = static_cast<size_t>(record);

        auto put_scalar = [&](const int varid, const double value)
        {
            const std::array<size_t, 1> start = {rec};
            const std::array<size_t, 1> count = {1};
            check(nc_put_vara_double(ncid, varid, start.data(), count.data(), &value), "scalar");
        };

        put_scalar(varid_time, time);
        put_scalar(varid_T, T);
        put_scalar(varid_p, p);
        put_scalar(varid_rho, rho);
        put_scalar(varid_rv, rv);
        put_scalar(varid_rs, rs);
        put_scalar(varid_supersaturation, rv/rs - 1.);

        std::array<double, Salsa::Gas::count> gas;
        for (int g=0; g<Salsa::Gas::count; ++g)
            gas[g] = fields(layout.gas_field(g), 0, 0, 0);

        const std::array<size_t, 2> start_gas = {rec, 0};
        const std::array<size_t, 2> count_gas = {1, Salsa::Gas::count};
        check(nc_put_vara_double(ncid, varid_gas, start_gas.data(), count_gas.data(), gas.data()), "gas");

        for (int c=0; c<3; ++c)
        {
            const size_t n_bins = category_end[c] - category_start[c];
            const size_t n_species = layout.n_species;

            const size_t n_dry = n_species - 1;
            const double number_limit = (c == 2) ? Salsa::prlim : Salsa::nlim;

            Kokkos::View<double*, Kokkos::HostSpace> number("number", n_bins);
            Kokkos::View<double**, Kokkos::HostSpace> mass("mass", n_species, n_bins);
            Kokkos::View<double*, Kokkos::HostSpace> diameter_wet("diameter_wet", n_bins);
            Kokkos::View<double*, Kokkos::HostSpace> diameter_dry("diameter_dry", n_bins);

            for (size_t b=0; b<n_bins; ++b)
            {
                const int bin = category_start[c] + b;
                number(b) = fields(layout.number_field(bin), 0, 0, 0);
                for (size_t s=0; s<n_species; ++s)
                    mass(s, b) = fields(layout.mass_field(bin, s), 0, 0, 0);

                // Volume-equivalent diameter of the mean particle, as the Fortran driver.
                if (number(b) > number_limit)
                {
                    double volume_wet = 0.;
                    double volume_dry = 0.;
                    for (size_t s=0; s<n_species; ++s)
                    {
                        volume_wet += mass(s, b)/layout.density[s];
                        if (s < n_dry)
                            volume_dry += mass(s, b)/layout.density[s];
                    }

                    diameter_wet(b) = std::pow(volume_wet/number(b)/Salsa::pi6, 1./3.);
                    diameter_dry(b) = std::pow(volume_dry/number(b)/Salsa::pi6, 1./3.);
                }
                else
                {
                    diameter_wet(b) = NC_FILL_DOUBLE;
                    diameter_dry(b) = NC_FILL_DOUBLE;
                }
            }

            const std::array<size_t, 2> start_bins = {rec, 0};
            const std::array<size_t, 2> count_bins = {1, n_bins};
            const std::array<size_t, 3> start_mass = {rec, 0, 0};
            const std::array<size_t, 3> count_mass = {1, n_species, n_bins};

            check(nc_put_vara_double(ncid, varid_number[c], start_bins.data(), count_bins.data(), number.data()), "number");
            check(nc_put_vara_double(ncid, varid_mass[c], start_mass.data(), count_mass.data(), mass.data()), "mass");
            check(nc_put_vara_double(ncid, varid_diameter_wet[c], start_bins.data(), count_bins.data(), diameter_wet.data()),
                    "diameter_wet");
            check(nc_put_vara_double(ncid, varid_diameter_dry[c], start_bins.data(), count_bins.data(), diameter_dry.data()),
                    "diameter_dry");
        }

        ++record;
    }
}
