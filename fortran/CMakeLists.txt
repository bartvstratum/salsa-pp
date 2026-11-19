# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Fortran reference build: vendored UCLALES-SALSA + new stand-ins and driver (fortran/driver).
# Note: the reference program links fortran/uclales-gpl (GPL v3+) and is therefore GPL v3 as a whole.
cmake_minimum_required(VERSION 3.20)
project(salsa_fortran_reference LANGUAGES Fortran)

if (NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "The Fortran reference is built with gfortran, as upstream.")
endif()

# Code generation flags as in upstream UCLALES-SALSA Makefile_ubuntu (double precision via -fdefault-real-8).
# The build type is ignored, so these are the only flags.
set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_RELEASE "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(upstream_flags -O2 -fdefault-real-8 -std=f2008 -fbounds-check -g3 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)

# Upstream also uses -Wall -Wtabs; warnings do not change the generated code, so they are only enabled for new code.
set(warning_flags -Wall -Wtabs)

file(GLOB vendored_sources CONFIGURE_DEPENDS
    uclales-salsa/src/*/*.f90
    uclales-gpl/src/*/*.f90)

set(stand_in_sources
    driver/grid.f90
    driver/util.f90)

add_library(uclales_salsa STATIC ${vendored_sources} ${stand_in_sources})
target_compile_options(uclales_salsa PRIVATE ${upstream_flags})
set_source_files_properties(driver/grid.f90 PROPERTIES COMPILE_OPTIONS "${warning_flags}")
# The stubs in util.f90 must keep upstream's interfaces but do not use their arguments.
set_source_files_properties(driver/util.f90 PROPERTIES COMPILE_OPTIONS "${warning_flags};-Wno-unused-dummy-argument")
set_target_properties(uclales_salsa PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
target_include_directories(uclales_salsa PUBLIC ${CMAKE_BINARY_DIR}/modules)

# SALSA++ git commit at configure time, written into the output files. Untracked files do not count as dirty.
find_package(Git QUIET)
set(GIT_COMMIT "unknown")
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if (NOT "${GIT_STATUS}" STREQUAL "")
        set(GIT_COMMIT "${GIT_COMMIT}-dirty")
    endif()
endif()
configure_file(driver/salsa_0d_version.f90.in ${CMAKE_BINARY_DIR}/generated/salsa_0d_version.f90 @ONLY)

# NetCDF-Fortran for the driver's input and output.
execute_process(COMMAND nf-config --includedir OUTPUT_VARIABLE NETCDF_FORTRAN_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
find_library(NETCDF_FORTRAN_LIB netcdff REQUIRED)
find_library(NETCDF_C_LIB netcdf REQUIRED)

add_executable(salsa_0d
    driver/salsa_0d.f90
    driver/salsa_0d_io.f90
    ${CMAKE_BINARY_DIR}/generated/salsa_0d_version.f90)
target_compile_options(salsa_0d PRIVATE ${upstream_flags} ${warning_flags})
target_include_directories(salsa_0d PRIVATE ${NETCDF_FORTRAN_INCLUDE})
target_link_libraries(salsa_0d PRIVATE uclales_salsa ${NETCDF_FORTRAN_LIB} ${NETCDF_C_LIB})
