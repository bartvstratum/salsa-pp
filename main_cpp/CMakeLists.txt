# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Drivers.

# SALSA++ git commit at configure time, written into the output files. Untracked files do not count as dirty.
find_package(Git QUIET)
set(GIT_COMMIT "unknown")
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if (NOT "${GIT_STATUS}" STREQUAL "")
        set(GIT_COMMIT "${GIT_COMMIT}-dirty")
    endif()
endif()
configure_file(salsa_version.h.in ${CMAKE_BINARY_DIR}/generated/salsa_version.h @ONLY)

add_executable(salsa_test salsa_test.cpp)
target_link_libraries(salsa_test PRIVATE salsa)

add_executable(salsa_0d salsa_0d.cpp salsa_0d_io.cpp)
target_include_directories(salsa_0d PRIVATE ${CMAKE_BINARY_DIR}/generated SYSTEM ${INCLUDE_DIRS})
target_link_libraries(salsa_0d PRIVATE salsa ${LIBS})
