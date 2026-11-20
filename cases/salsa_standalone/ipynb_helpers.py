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


def load_runs(case):
    case_dir = pathlib.Path(__file__).resolve().parent / f'case_{case}'
    runs = {'Fortran': xr.open_dataset(case_dir / 'output_fortran.nc')}
    for path in sorted(case_dir.glob('output_cpp*.nc')):
        runs[path.stem.removeprefix('output_')] = xr.open_dataset(path)
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


def title(ds, c):
    b = regimes(ds, c).get('b')
    used = b is not None and float(ds[f'number_{c}'].isel({f'bin_{c}': b}).max()) > 0
    return c + (' (light: regime b)' if used else '')


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
                d = ds[f'diameter_wet_{c}'].sel(time=t).values * 1e6
                n = ds[f'number_{c}'].sel(time=t).values * 1e-6
                for regime, bins in regimes(ds, c).items():
                    label = f't = {t} s' if name == 'Fortran' and regime != 'b' else None
                    ax.plot(d[bins], n[bins], line_style(name), marker='o', color=color,
                            alpha=0.5 if regime == 'b' else 1, label=label)
        ax.set_xscale('log')
        ax.set_title(title(runs['Fortran'], c))
        ax.grid()
    axs[0,0].legend(ncol=len(times))
    axs[-1,0].set_xlabel('wet diameter (µm)')
    fig.supylabel('number per bin (mg$^{-1}$)')
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
                    ax.plot(ds.time, ds[f'diameter_wet_{c}'].isel({f'bin_{c}': b}) * 1e6, line_style(name),
                            color=color, alpha=0.5 if regime == 'b' else 1)
        ax.set_yscale('log')
        ax.set_title(title(ref, c))
        ax.grid()
    axs[-1,0].set_xlabel('time (s)')
    fig.supylabel('wet diameter (µm)')
    fig.tight_layout()
    plt.show()


def plot_conservation(runs):
    fig, ax = plt.subplots()
    for name, ds in runs.items():
        sp = species(ds)
        totals = {'number': sum(ds[f'number_{c}'].sum(f'bin_{c}') for c in categories),
                  'water (vapour + liquid)': ds.qv + sum(liquid_water(ds, c) for c in categories)}
        for i, s in enumerate(sp):
            if s != 'H2O':
                totals[f'{s} mass'] = sum(ds[f'mass_{c}'].isel(species=i).sum(f'bin_{c}') for c in categories)
        for label, x in totals.items():
            if x[0] > 0:
                ax.plot(ds.time, x / x[0] - 1, line_style(name), label=f'{label}, {name}')
    ax.set_xlabel('time (s)')
    ax.set_ylabel('relative change since t = 0 (-)')
    ax.legend()
    ax.grid()
    fig.tight_layout()
    plt.show()


def compare(runs):
    ref = runs['Fortran']
    rows = []
    for name, ds in runs.items():
        if name == 'Fortran':
            continue
        for var in ref.data_vars:
            if 'time' not in ref[var].dims:
                continue
            diff = float(abs(ds[var] - ref[var]).max())
            scale = float(abs(ref[var]).max())
            rows.append((name, var, diff, diff / scale if scale > 0 else np.nan))
    if not rows:
        print('No C++ output yet.')
        return None
    return pd.DataFrame(rows, columns=['run', 'variable', 'max |diff|', 'max |diff| / max |Fortran|'])
