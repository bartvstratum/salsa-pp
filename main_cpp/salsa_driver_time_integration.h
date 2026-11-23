/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Translated from UCLALES-SALSA (commit 8c77958), constrain_SALSA.f90.
 * Copyright (c) 2020 Finnish Meteorological Institute, MIT licence (see LICENSE-MIT).
 *
 * Time integration of the standalone drivers, as the LES side of UCLALES-SALSA: limit the SALSA tendencies and
 * update the fields. Not part of the solver; a host model does this itself.
 */

#pragma once

#include "salsa_types.h"


namespace Salsa
{
    // Fortran: constrain_SALSA.f90: tend_constrain2
    // A step lowers a field to at most 1e-8 of its old value.
    inline void constrain_tendencies(
            const Array_4d<const TF>& fields,
            const Array_4d<TF>& fields_tend,
            const TF dt,
            const int istart,
            const int iend,
            const int jstart,
            const int jend,
            const int kstart,
            const int kend)
    {
        const int n_fields = fields.extent(0);

        Kokkos::parallel_for("constrain_tendencies",
                Range_3d({kstart, jstart, istart}, {kend, jend, iend}),
                KOKKOS_LAMBDA(const int k, const int j, const int i)
                {
                    for (int n=0; n<n_fields; ++n)
                        fields_tend(n, k, j, i) = Kokkos::max((-1._tf + 1.e-8_tf)*fields(n, k, j, i)/dt, fields_tend(n, k, j, i));
                });
    }


    // Forward Euler update of the LES, as the Fortran driver (salsa_0d.f90). The initialisation call clips the
    // fields at zero (init.f90). The tendencies are reset afterwards.
    inline void integrate_forward_euler(
            const Array_4d<TF>& fields,
            const Array_4d<TF>& fields_tend,
            const Array_3d<TF>& rv,
            const Array_3d<TF>& rv_tend,
            const TF dt,
            const bool clip_at_zero,
            const int istart,
            const int iend,
            const int jstart,
            const int jend,
            const int kstart,
            const int kend)
    {
        const int n_fields = fields.extent(0);

        Kokkos::parallel_for("integrate_forward_euler",
                Range_3d({kstart, jstart, istart}, {kend, jend, iend}),
                KOKKOS_LAMBDA(const int k, const int j, const int i)
                {
                    for (int n=0; n<n_fields; ++n)
                    {
                        if (clip_at_zero)
                            fields(n, k, j, i) = Kokkos::max(fields(n, k, j, i) + dt*fields_tend(n, k, j, i), 0._tf);
                        else
                            fields(n, k, j, i) = fields(n, k, j, i) + dt*fields_tend(n, k, j, i);

                        fields_tend(n, k, j, i) = 0._tf;
                    }

                    rv(k, j, i) += dt*rv_tend(k, j, i);
                    rv_tend(k, j, i) = 0._tf;
                });
    }
}
