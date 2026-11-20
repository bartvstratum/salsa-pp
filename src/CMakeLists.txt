# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# The solver library, without I/O.
add_library(salsa STATIC
    salsa_layout.cpp
    salsa_solver.cpp)

target_include_directories(salsa PUBLIC "../include")
target_link_libraries(salsa PUBLIC Kokkos::kokkos)
