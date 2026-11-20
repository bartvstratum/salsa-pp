# SALSA++
C++/Kokkos port of the sectional aerosol–cloud–precipitation microphysics SALSA (Sectional Aerosol module for Large Scale Applications), as extended in UCLALES-SALSA.

Licence: Apache 2.0, except the vendored third-party code in `fortran/`, which keeps its own licence (see the README in each directory there).

## SALSA++
Requires a C++20 compiler and CMake ≥ 3.20. Kokkos is a submodule (`git submodule update --init`). System settings are in `config/<system>.cmake`, as in MicroHH:

    mkdir -p build/cpu
    cd build/cpu
    cmake ../.. -DSYST=ubuntu_22lts_gcc
    make -j
    ./main_cpp/salsa_test

One build directory per configuration, e.g. `build/gpu_sp` with `cmake ../.. -DSYST=ubuntu_22lts_gcc -DUSEGPU=1 -DUSESP=1`.

Options: `-DUSESP=1` for single precision, `-DUSEGPU=1` for GPU. With CUDA, Kokkos needs nvcc ≥ 12.2, first in `PATH` or via `CUDA_ROOT` (e.g. `CUDA_ROOT=/usr/local/cuda`).

## Fortran reference
Requires gfortran, CMake and NetCDF-Fortran:

    mkdir -p build/fortran
    cd build/fortran
    cmake ../../fortran
    make -j

## Reference cases
The case inputs are created with Python:

    python -m venv /path/to/venv
    source /path/to/venv/bin/activate
    pip install -r requirements.txt
    python cases/salsa_standalone/salsa_standalone_input.py

This writes `input.nc` in `cases/salsa_standalone/case_a` to `case_d`. Run a case from its directory, which holds the `NAMELIST`:

    cd cases/salsa_standalone/case_a
    ../../../build/fortran/salsa_0d
    ../../../build/cpu/main_cpp/salsa_0d

The Fortran output goes to `output_fortran.nc`, the C++ output to `output_cpp.nc`. The notebooks `cases/salsa_standalone/compare_case_*.ipynb` plot the cases.
