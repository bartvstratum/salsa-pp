# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# The solver library, without I/O.
add_library(salsa STATIC
    salsa_layout.cpp
    salsa_sizedist.cpp
    salsa_solver.cpp)

target_include_directories(salsa PUBLIC "../include" "../include_kernels")
target_link_libraries(salsa PUBLIC Kokkos::kokkos)
