# SALSA++
C++/Kokkos port of the sectional aerosol–cloud–precipitation microphysics SALSA (Sectional Aerosol module for Large Scale Applications), as extended in UCLALES-SALSA.

Licence: Apache 2.0, except the vendored third-party code in `fortran/`, which keeps its own licence (see the README in each directory there).

## Fortran to C++/Kokkos porting
Most of the translation was done by Claude Opus 5.5 under human supervision.

- The numerics and physics closely follow the original Fortran code, but the code design differs substantially: it targets performance on CPUs and GPUs with support for single precision.
- First, all processes were ported one-to-one, deliberately leaving existing bugs and single-precision issues of the original code in place.
- When all processes up to test case D had been ported, this state was tagged as [v0.1.0](https://github.com/bartvstratum/salsa-pp/releases/tag/v0.1.0). Up to that point, double-precision results from the Fortran and C++ code were exactly identical.
- After v0.1.0, issues were fixed step by step in both codes, keeping the Fortran/C++ tests reproducible. Each fix has its own GitHub [issue](https://github.com/bartvstratum/salsa-pp/issues), so the diff for each change is well documented.

## SALSA++
Requires a C++20 compiler and CMake ≥ 3.20. Kokkos is a submodule (`git submodule update --init`). System settings are in `config/<system>.cmake`, as in MicroHH:

    mkdir -p build/cpu
    cd build/cpu
    cmake ../.. -DSYST=ubuntu_22lts_gcc
    make -j
    ./main_cpp/salsa_bench 32 32 32 1

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
