# SALSA++
C++/Kokkos port of the sectional aerosol–cloud–precipitation microphysics SALSA (Sectional Aerosol module for Large Scale Applications), as extended in UCLALES-SALSA.

Licence: Apache 2.0, except the vendored third-party code in `fortran/`, which keeps its own licence (see the README in each directory there).

## Fortran reference
Requires gfortran, CMake and NetCDF-Fortran. From the repository root:

    cmake -S fortran -B build/fortran
    cmake --build build/fortran -j

## Reference cases
The case inputs are created with Python:

    python -m venv /path/to/venv
    source /path/to/venv/bin/activate
    pip install -r requirements.txt
    python cases/salsa_standalone/salsa_standalone_input.py

This writes `input.nc` in `cases/salsa_standalone/case_a`, `case_b` and `case_c`. Run a case from its directory, which holds the `NAMELIST`:

    cd cases/salsa_standalone/case_a
    ../../../build/fortran/salsa_0d

The output goes to `output_fortran.nc`.
