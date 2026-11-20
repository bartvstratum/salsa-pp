/*
 * SALSA++
 * Copyright (c) 2026 Wageningen University & Research
 * SPDX-License-Identifier: Apache-2.0
 *
 * Floating point type, execution space and array types, as in MicroHH.
 */

#pragma once

#include <Kokkos_Core.hpp>


namespace Salsa
{
    #ifdef USE_SINGLE_PRECISION
    using TF = float;
    #else
    using TF = double;
    #endif

    KOKKOS_INLINE_FUNCTION constexpr TF operator""_tf(const long double ld) { return static_cast<TF>(ld); }
    KOKKOS_INLINE_FUNCTION constexpr TF operator""_tf(const unsigned long long ll) { return static_cast<TF>(ll); }


    // Default execution space: Cuda/HIP for GPU builds, OpenMP or Serial for CPU.
    #if defined(USECUDA)
    using Default_exec = Kokkos::Cuda;
    #elif defined(USEHIP)
    using Default_exec = Kokkos::HIP;
    #elif defined(KOKKOS_ENABLE_OPENMP)
    using Default_exec = Kokkos::OpenMP;
    #else
    using Default_exec = Kokkos::Serial;
    #endif

    using Range_3d = Kokkos::MDRangePolicy<Default_exec, Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;


    // C-style index order, no overlap in memory.
    template<typename T>
    using Array_3d = Kokkos::View<T***, Kokkos::LayoutRight, Default_exec::memory_space, Kokkos::MemoryTraits<Kokkos::Restrict>>;

    template<typename T>
    using Array_4d = Kokkos::View<T****, Kokkos::LayoutRight, Default_exec::memory_space, Kokkos::MemoryTraits<Kokkos::Restrict>>;
}
