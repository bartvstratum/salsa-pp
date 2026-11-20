# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Settings for Ubuntu with GCC, and CUDA for GPU builds. Kokkos needs nvcc >= 12.2: put it first in PATH,
# or set CUDA_ROOT (e.g. CUDA_ROOT=/usr/local/cuda), which Kokkos' nvcc_wrapper uses.

# Set the compiler via environment variables.
if (USEGPU)
    set(ENV{CXX} ${CMAKE_CURRENT_LIST_DIR}/../extern/kokkos/bin/nvcc_wrapper)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
else()
    set(ENV{CXX} g++)
endif()

# No FMA contraction (-ffp-contract=off): with -march=native, g++ fuses a*b + c into one instruction with different
# rounding than the Fortran reference, which is built without -march. Without FMA the C++ results are bit-identical.
set(USER_CXX_FLAGS "-Wall -fPIC -fopenmp -ffp-contract=off")
set(USER_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
set(USER_CXX_FLAGS_DEBUG "-g -O0")

# Not forced, so that -DKokkos_ENABLE_OPENMP=OFF gives a Serial build.
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)


# Set the libraries.
set(NETCDF_LIB_C "netcdf")
set(LIBS ${NETCDF_LIB_C})
