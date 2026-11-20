# SALSA++
# Copyright (c) 2026 Wageningen University & Research
# SPDX-License-Identifier: Apache-2.0
#
# Plots for the compare_case_*.ipynb notebooks.

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

categories = ['aerosol', 'cloud', 'rain']

# Identical results hide each other, so the size distribution and wet diameter plots multiply the C++/Kokkos diameters
# by this factor.
diameter_shift = 1.1


def load_runs(case):
    case_dir = pathlib.Path(__file__).resolve().parent / f'case_{case}'
    runs = {'Fortran': xr.open_dataset(case_dir / 'output_fortran.nc')}
    for path in sorted(case_dir.glob('output_cpp*.nc')):
        name = path.stem.removeprefix('output_').replace('cpp', 'C++/Kokkos').replace('_', ' ')
        runs[name] = xr.open_dataset(path)
    return runs


def species(ds):
    return [s.decode().strip() if isinstance(s, bytes) else str(s).strip() for s in ds['species'].values]


def non_empty(runs):
    return [c for c in categories if any(float(ds[f'number_{c}'].max()) > 0 for ds in runs.values())]


def line_style(name):
    return '-' if name == 'Fortran' else '--'


# Regime b bins follow the regime a bins; the bin centre diameter drops where b starts.
def regimes(ds, c):
    d = ds[f'bin_{c}'].values
    reset = np.where(np.diff(d) < 0)[0]
    if reset.size == 0:
        return {'': slice(None)}
    return {'a': slice(0, reset[0] + 1), 'b': slice(reset[0] + 1, None)}


# One panel per row with a shared x-axis; width from matplotlibrc, taller with more panels.
def stacked(n, **kwargs):
    w, h = plt.rcParams['figure.figsize']
    return plt.subplots(n, 1, squeeze=False, sharex=True, figsize=(w, h * (1 + 0.2 * (n - 1))), **kwargs)


def diameter_factor(name):
    return 1 if name == 'Fortran' else diameter_shift


def title(runs, c):
    ds = runs['Fortran']
    b = regimes(ds, c).get('b')
    notes = []
    if b is not None and float(ds[f'number_{c}'].isel({f'bin_{c}': b}).max()) > 0:
        notes.append('light: regime b')
    if len(runs) > 1:
        notes.append(f'C++/Kokkos dashed, diameter shifted by +{(diameter_shift - 1) * 100:.0f} %')
    if not notes:
        return c
    return c + ' (' + '; '.join(notes) + ')'


def bin_range(ds, c, b):
    low, high = ds[f'bin_{c}_bounds'].values[b] * 1e6
    return f'{low:.3g}-{high:.3g} µm'


def liquid_water(ds, c):
    return ds[f'mass_{c}'].isel(species=species(ds).index('H2O')).sum(f'bin_{c}')


def summary(runs):
    ref = runs['Fortran']
    print(ref.attrs['title'])
    print('Runs:     ', list(runs))
    print('Species:  ', species(ref))
    print('Bins:     ', {c: ref.sizes[f'bin_{c}'] for c in categories})
    print('Not empty:', non_empty(runs))


def plot_thermodynamics(runs):
    fig, axs = plt.subplots(2, 2, sharex=True)
    for name, ds in runs.items():
        ls = line_style(name)
        axs[0,0].plot(ds.time, ds.temperature, ls, label=name)
        axs[0,1].plot(ds.time, ds.supersaturation * 100, ls, label=name)
        axs[1,0].plot(ds.time, ds.qv * 1e3, ls, label=name)
        for c in non_empty(runs):
            axs[1,1].plot(ds.time, liquid_water(ds, c) * 1e3, ls, label=f'{c}, {name}')
    axs[0,0].set_ylabel('temperature (K)')
    axs[0,1].set_ylabel('supersaturation (%)')
    axs[1,0].set_ylabel('water vapour (g kg$^{-1}$)')
    axs[1,1].set_ylabel('liquid water, sum over bins (g kg$^{-1}$)')
    for ax in axs.flat:
        ax.legend()
        ax.grid()
    for ax in axs[1]:
        ax.set_xlabel('time (s)')
    fig.tight_layout()
    plt.show()


def plot_size_distributions(runs, times=(0, 100, 200, 300, 400)):
    cats = non_empty(runs)
    colors = [f'C{i}' for i in range(len(times))]
    fig, axs = stacked(len(cats))
    for ax, c in zip(axs[:,0], cats):
        for name, ds in runs.items():
            for t, color in zip(times, colors):
                if t not in ds.time:
                    continue
                d = ds[f'diameter_wet_{c}'].sel(time=t).values * 1e6 * diameter_factor(name)
                n = ds[f'number_{c}'].sel(time=t).values * 1e-6
                for regime, bins in regimes(ds, c).items():
                    label = f't = {t} s' if name == 'Fortran' and regime != 'b' else None
                    ax.plot(d[bins], n[bins], line_style(name), marker='o', color=color,
                            alpha=0.5 if regime == 'b' else 1, label=label)
        ax.set_xscale('log')
        ax.set_title(title(runs, c))
        ax.set_ylabel('number per bin (mg$^{-1}$)')
        ax.grid()
    axs[0,0].legend(ncol=min(len(times), 5))
    axs[-1,0].set_xlabel('wet diameter (µm)')
    fig.tight_layout()
    plt.show()


def plot_wet_diameters(runs):
    cats = non_empty(runs)
    fig, axs = stacked(len(cats), sharey=True)
    for ax, c in zip(axs[:,0], cats):
        ref = runs['Fortran']
        for regime, bins in regimes(ref, c).items():
            idx = np.arange(ref.sizes[f'bin_{c}'])[bins]
            colors = plt.cm.plasma(np.linspace(0, 0.9, idx.size))
            for name, ds in runs.items():
                for b, color in zip(idx, colors):
                    d = ds[f'diameter_wet_{c}'].isel({f'bin_{c}': b}).values * 1e6 * diameter_factor(name)
                    alpha = 0.5 if regime == 'b' else 1
                    ax.plot(ds.time, d, line_style(name), color=color, alpha=alpha)

                    # Dry size range of the bin at the end of the line.
                    valid = np.where(np.isfinite(d))[0]
                    if name == 'Fortran' and valid.size > 0:
                        last = valid[-1]
                        ax.text(ds.time.values[last], d[last], ' ' + bin_range(ds, c, b), color=color, alpha=alpha,
                                fontsize='x-small', va='center')
        ax.set_yscale('log')
        ax.set_title(title(runs, c))
        ax.set_ylabel('wet diameter (µm)')
        ax.grid()
    t = runs['Fortran'].time.values
    axs[-1,0].set_xlim(t[0], t[-1] + 0.25 * (t[-1] - t[0]))
    axs[-1,0].set_xlabel('time (s)')
    fig.tight_layout()
    plt.show()


# Water and the dry compounds are conserved; the total number is not, since collisions merge particles.
def plot_conservation(runs):
    fig, axs = stacked(2)
    for name, ds in runs.items():
        ls = line_style(name)
        totals = {'water (vapour + liquid)': ds.qv + sum(liquid_water(ds, c) for c in categories)}
        for i, s in enumerate(species(ds)):
            if s != 'H2O':
                totals[s] = sum(ds[f'mass_{c}'].isel(species=i).sum(f'bin_{c}') for c in categories)
        for label, x in totals.items():
            if x[0] > 0:
                axs[0,0].plot(ds.time, x / x[0] - 1, ls, label=f'{label}, {name}')
        n = sum(ds[f'number_{c}'].sum(f'bin_{c}') for c in categories)
        axs[1,0].plot(ds.time, n / n[0] - 1, ls, label=name)
    axs[0,0].set_title('mass, conserved')
    axs[1,0].set_title('number, decreases by collisions')
    for ax in axs[:,0]:
        ax.set_ylabel('relative change since t = 0 (-)')
        ax.legend()
        ax.grid()
    axs[-1,0].set_xlabel('time (s)')
    fig.tight_layout()
    plt.show()


# Compares by position, not by coordinate value. A bin that is empty (NaN) in one run only counts as an infinite difference.
def compare(runs):
    ref = runs['Fortran']
    rows = []
    for name, ds in runs.items():
        if name == 'Fortran':
            continue
        # All variables, as xarray makes gas (named as its dimension) a coordinate.
        for var in ref.variables:
            if var == 'time' or 'time' not in ref[var].dims:
                continue
            a = ref[var].values
            b = ds[var].values
            diff = np.where(np.isnan(a) != np.isnan(b), np.inf, np.abs(b - a))
            diff = float(np.nanmax(diff)) if not np.all(np.isnan(diff)) else 0.
            scale = float(np.nanmax(np.abs(a))) if not np.all(np.isnan(a)) else 0.
            rows.append((name, var, diff, diff / scale if scale > 0 else np.nan))
    if not rows:
        print('No C++ output yet.')
        return None
    return pd.DataFrame(rows, columns=['run', 'variable', 'max |diff|', 'max |diff| / max |Fortran|'])
