# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Drivers.
add_executable(salsa_test salsa_test.cpp)
target_link_libraries(salsa_test PRIVATE salsa)
