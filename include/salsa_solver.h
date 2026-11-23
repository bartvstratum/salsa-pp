/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * The solver: one SALSA step for all grid points.
 */

#pragma once

#include <array>
#include <string>
#include <vector>

#include "salsa_types.h"
#include "salsa_layout.h"


namespace Salsa
{
    class Solver
    {
        public:
            Solver(
                    const std::array<TF, 4>& regime_limits,         // Fortran: reglim
                    const std::array<int, 2>& n_bins_regime,        // Fortran: nbin
                    const int n_bins_rain,                          // Fortran: bloPrc%nbins
                    const TF rain_diameter_low,                     // Fortran: bloPrc%dlo
                    const TF rain_volume_ratio,                     // Fortran: bloPrc%vol_ratio
                    const std::vector<std::string>& dry_species,    // Fortran: listspec
                    const int coagulation_mode,                     // Fortran: lscoag%mode
                    const int autoconversion_mode,                  // Fortran: lsauto%mode
                    const int activation_mode,                      // Fortran: lsactiv%mode
                    const bool coagulate_aer,                       // Fortran: lscgaa
                    const bool coagulate_cld,                       // Fortran: lscgcc
                    const bool cld_collects_aer,                    // Fortran: lscgca
                    const bool rain_collects_cld,                   // Fortran: lscgpc
                    const bool rain_collects_aer,                   // Fortran: lscgpa
                    const bool coagulate_rain,                      // Fortran: lscgpp
                    const bool condense_gases,                      // Fortran: lscndgas
                    const bool condense_water_aer,                  // Fortran: lscndh2oae
                    const bool condense_water_cld,                  // Fortran: lscndh2ocl
                    const TF rh_limit);                             // Fortran: rhlim

            Solver(const Solver&) = delete;
            Solver& operator=(const Solver&) = delete;
            Solver(Solver&&) = delete;
            Solver& operator=(Solver&&) = delete;

            const Layout& get_layout() const
            {
                return layout;
            }

            // Runs on the grid points [istart, iend) x [jstart, jend) x [kstart, kend), and adds the change of the
            // fields and of rv over dt as tendencies, as run_SALSA in UCLALES-SALSA. The host applies them.
            void exec(
                    const Array_4d<TF>& fields_tend,                // (n, k, j, i), per s
                    const Array_3d<TF>& rv_tend,                    // kg kg-1 s-1
                    const Array_4d<TF>& diameter_wet,               // (bin, k, j, i), kept between calls as in UCLALES-SALSA
                    const Array_4d<const TF>& fields,               // (n, k, j, i), n in the order of the layout
                    const Array_3d<const TF>& rv,                   // Water vapour mixing ratio (kg kg-1)
                    const Array_3d<const TF>& T,                    // K
                    const Array_3d<const TF>& p,                    // Pa
                    const Array_3d<const TF>& rs,                   // Saturation mixing ratio over water (kg kg-1)
                    const Array_3d<const TF>& rho,                  // kg m-3
                    const TF dt,
                    const bool initialize,                          // Equilibrate the aerosol water first; Fortran: run_SALSA initialize
                    const bool do_coagulation,
                    const bool update_coagulation_kernels,          // Coagulation mode 2 only
                    const bool do_condensation,
                    const bool do_activation,
                    const bool do_autoconversion,
                    const bool do_update_distribution,              // Fortran: lsdistupdate
                    const bool free_rh,                             // Otherwise RH is limited to rh_limit
                    const int istart,
                    const int iend,
                    const int jstart,
                    const int jend,
                    const int kstart,
                    const int kend);

            // Fortran: constrain_SALSA.f90: SALSA_diagnostics (without clipping)
            // Removes small particles on the grid points [istart, iend) x [jstart, jend) x [kstart, kend), on the
            // fields after the time integration, as UCLALES-SALSA after each step.
            void remove_small_particles(
                    const Array_4d<TF>& fields,                     // (n, k, j, i), n in the order of the layout
                    const Array_3d<const TF>& T,                    // K
                    const Array_3d<const TF>& rh,                   // Relative humidity at the start of the step
                    const int istart,
                    const int iend,
                    const int jstart,
                    const int jend,
                    const int kstart,
                    const int kend);

        private:
            const Layout layout;

            Array_1d<TF> volume_low;
            Array_1d<TF> volume_high;
            Array_1d<TF> density;
            Array_1d<TF> molar_mass;
            Array_1d<TF> diameter_mid;
            Array_1d<TF> number_limit;
            Array_1d<TF> diameter_limit;
            Array_1d<TF> diameter_low;                      // Diameter of the lower volume limit; Fortran: precpbins (rain)
            Array_1d<TF> dissociation;
            Array_1d<TF> volume_ratio_low;
            Array_1d<TF> volume_ratio_high;
            Array_1d<bool> soluble;

            const int coagulation_mode;
            const int autoconversion_mode;
            const int activation_mode;

            const bool coagulate_aer;
            const bool coagulate_cld;
            const bool cld_collects_aer;
            const bool rain_collects_cld;
            const bool rain_collects_aer;
            const bool coagulate_rain;

            const bool condense_gases;
            const bool condense_water_aer;
            const bool condense_water_cld;

            const TF rh_limit;

            // Species indices, -1 if not used.
            const int iso4;
            const int ioc;
            const int ibc;
            const int idu;
            const int iss;
            const int iwa;
            const int ino;
            const int inh;
    };
}
