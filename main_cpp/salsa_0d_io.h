/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * NetCDF input and output of the 0D driver. The output has the layout of the Fortran driver's output.
 */

#pragma once

#include <array>
#include <string>
#include <vector>

#include "salsa_types.h"
#include "salsa_layout.h"


namespace Salsa_0d
{
    using Array_4d_h = Salsa::Array_4d<Salsa::TF>::host_mirror_type;


    std::vector<double> read_variable(const std::string& file_name, const std::string& name);
    std::vector<std::string> read_names(const std::string& file_name, const std::string& name);


    class Output
    {
        public:
            Output(
                    const std::string& file_name,
                    const std::string& case_name,
                    const std::string& input_file,
                    const std::string& namelist_text,
                    const Salsa::Layout& layout);

            ~Output();

            Output(const Output&) = delete;
            Output& operator=(const Output&) = delete;
            Output(Output&&) = delete;
            Output& operator=(Output&&) = delete;

            // The fields of the grid point (0, 0, 0) as the next time record.
            void write(
                    const Array_4d_h& fields,
                    const double time,
                    const double T,
                    const double p,
                    const double rho,
                    const double rv,
                    const double rs);

        private:
            const Salsa::Layout& layout;

            // Aerosol, cloud and rain.
            std::array<int, 3> category_start;
            std::array<int, 3> category_end;

            int ncid;
            int record;

            int varid_time;
            int varid_T;
            int varid_p;
            int varid_rho;
            int varid_rv;
            int varid_rs;
            int varid_supersaturation;
            int varid_gas;

            std::array<int, 3> varid_number;
            std::array<int, 3> varid_mass;
            std::array<int, 3> varid_diameter_wet;
            std::array<int, 3> varid_diameter_dry;
    };
}
