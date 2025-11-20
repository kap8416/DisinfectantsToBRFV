#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Advanced Emax dose–response analysis for ToBRFV disinfectant efficacy

Features:
- Emax model with bounds (E0 and Emax constrained to [0, 100])
- Robust initial guesses based on data
- Bootstrap confidence intervals for E0, Emax, ED50 and ED90
- Monotonicity diagnostics for each surface–method combination
- Matplotlib figures (journal-ready) with 2x2 layout
- Export of parameter table to CSV (MDPI-friendly)

Author: Katia Aviña Padilla 
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------------------------------------------------
# 1. Configuration
# -------------------------------------------------------------------

DATA_FILE = "dose_response_cleaned.csv"  #
OUTPUT_TABLE = "emax_parameters_bootstrap.csv"
FIGURE_FILE = "Figure9_Emax_dose_response.png"

N_BOOT = 1000  # número de réplicas bootstrap (puedes bajar a 500 si tarda)

# Column names expected in DATA_FILE
COL_SURFACE = "surface"
COL_METHOD = "method"
COL_DOSE = "dose"
COL_EFF = "efficacy"

# -------------------------------------------------------------------
# 2. Emax model and helper functions
# -------------------------------------------------------------------

def emax_model(dose, E0, Emax, ED50):
    """
    Standard 3-parameter Emax model:

        E(d) = E0 + (Emax - E0) * dose / (ED50 + dose)

    E0   : response at dose = 0
    Emax : maximum achievable effect
    ED50 : dose producing 50% of the Emax–E0 effect

    We later derive ED90 analytically from ED50.
    """
    dose = np.asarray(dose)
    return E0 + (Emax - E0) * dose / (ED50 + dose)


def ed90_from_ed50(ED50):
    """
    For the 3-parameter Emax model, ED90 = 9 * ED50.
    This comes from setting effect to 90% of (Emax - E0).
    """
    return 9.0 * ED50


def initial_guesses(dose, eff):
    """
    Reasonable automatic initial guesses based on data.
    """
    E0_init = np.percentile(eff, 10)  # lower-end efficacy
    Emax_init = np.percentile(eff, 90)  # upper-end efficacy
    if Emax_init < E0_init:
        # fallback if data are noisy
        Emax_init = eff.max()
        E0_init = eff.min()

    # ED50: some dose in the mid-range (avoid zero)
    pos_dose = dose[dose > 0]
    if len(pos_dose) == 0:
        ED50_init = 1.0
    else:
        ED50_init = np.median(pos_dose)

    return E0_init, Emax_init, ED50_init


def fit_emax_group(dose, eff):
    """
    Fit the Emax model for a single group (surface–method).
    Returns:
        params  : dict with E0, Emax, ED50, ED90
        success : bool, whether the fit converged
    """
    # Convert to numpy arrays
    dose = np.asarray(dose, dtype=float)
    eff = np.asarray(eff, dtype=float)

    # Remove missing / NaN
    mask = ~np.isnan(dose) & ~np.isnan(eff)
    dose = dose[mask]
    eff = eff[mask]

    if len(dose) < 4:
        return None, False  # not enough data points

    # Initial guesses
    p0 = initial_guesses(dose, eff)

    # Bounds:
    # E0 in [0, 100], Emax in [0, 100], ED50 > 0 (here [1e-6, inf))
    lower_bounds = [0.0, 0.0, 1e-6]
    upper_bounds = [100.0, 100.0, np.inf]

    try:
        popt, pcov = curve_fit(
            emax_model,
            dose,
            eff,
            p0=p0,
            bounds=(lower_bounds, upper_bounds),
            maxfev=10000,
        )
        E0_hat, Emax_hat, ED50_hat = popt
        ED90_hat = ed90_from_ed50(ED50_hat)

        params = {
            "E0": E0_hat,
            "Emax": Emax_hat,
            "ED50": ED50_hat,
            "ED90": ED90_hat,
        }
        return params, True

    except Exception:
        return None, False


def bootstrap_emax(dose, eff, n_boot=N_BOOT, random_state=42):
    """
    Nonparametric bootstrap for Emax parameters.
    Returns dict with point estimates and 95% CIs.
    """
    rng = np.random.default_rng(random_state)

    params_hat, ok = fit_emax_group(dose, eff)
    if not ok:
        return None

    E0_list = []
    Emax_list = []
    ED50_list = []
    ED90_list = []

    n = len(dose)
    idx = np.arange(n)

    for _ in range(n_boot):
        # resample indices with replacement
        bs_idx = rng.choice(idx, size=n, replace=True)
        dose_bs = np.array(dose)[bs_idx]
        eff_bs = np.array(eff)[bs_idx]

        params_bs, ok_bs = fit_emax_group(dose_bs, eff_bs)
        if not ok_bs:
            continue

        E0_list.append(params_bs["E0"])
        Emax_list.append(params_bs["Emax"])
        ED50_list.append(params_bs["ED50"])
        ED90_list.append(params_bs["ED90"])

    def ci(series):
        if len(series) == 0:
            return np.nan, np.nan
        lower = np.percentile(series, 2.5)
        upper = np.percentile(series, 97.5)
        return lower, upper

    E0_L, E0_U = ci(E0_list)
    Emax_L, Emax_U = ci(Emax_list)
    ED50_L, ED50_U = ci(ED50_list)
    ED90_L, ED90_U = ci(ED90_list)

    out = {
        "E0": params_hat["E0"],
        "E0_L": E0_L,
        "E0_U": E0_U,
        "Emax": params_hat["Emax"],
        "Emax_L": Emax_L,
        "Emax_U": Emax_U,
        "ED50": params_hat["ED50"],
        "ED50_L": ED50_L,
        "ED50_U": ED50_U,
        "ED90": params_hat["ED90"],
        "ED90_L": ED90_L,
        "ED90_U": ED90_U,
    }
    return out


def check_monotonicity(dose, eff):
    """
    Simple monotonicity diagnostic:
    - Sort by dose
    - Compute Spearman correlation between dose and efficacy
    Returns correlation and a qualitative flag.
    """
    from scipy.stats import spearmanr

    dose = np.asarray(dose)
    eff = np.asarray(eff)
    mask = ~np.isnan(dose) & ~np.isnan(eff)
    dose = dose[mask]
    eff = eff[mask]

    if len(dose) < 4:
        return np.nan, "insufficient"

    r, p = spearmanr(dose, eff)
    if np.isnan(r):
        flag = "indeterminate"
    elif r > 0.3:
        flag = "monotone-increasing"
    elif r < -0.3:
        flag = "monotone-decreasing"
    else:
        flag = "weak/none"

    return r, flag


# -------------------------------------------------------------------
# 3. Load data
# -------------------------------------------------------------------

df = pd.read_csv(DATA_FILE)

# Ensure types
df[COL_DOSE] = pd.to_numeric(df[COL_DOSE], errors="coerce")
df[COL_EFF] = pd.to_numeric(df[COL_EFF], errors="coerce")

# Unique combinations
groups = df[[COL_SURFACE, COL_METHOD]].drop_duplicates()

results = []

# -------------------------------------------------------------------
# 4. Fit Emax model + bootstrap for each surface–method
# -------------------------------------------------------------------

for _, row in groups.iterrows():
    surface = row[COL_SURFACE]
    method = row[COL_METHOD]

    sub = df[(df[COL_SURFACE] == surface) & (df[COL_METHOD] == method)]
    dose = sub[COL_DOSE].values
    eff = sub[COL_EFF].values

    print(f"Fitting Emax for {surface} – {method} (n={len(sub)})")

    mono_r, mono_flag = check_monotonicity(dose, eff)

    boot_params = bootstrap_emax(dose, eff, n_boot=N_BOOT)
    if boot_params is None:
        print(f"  WARNING: Emax fit failed for {surface} – {method}")
        continue

    res = {
        "surface": surface,
        "method": method,
        "monotonicity_r": mono_r,
        "monotonicity_flag": mono_flag,
    }
    res.update(boot_params)
    results.append(res)

# Convert to DataFrame and save
res_df = pd.DataFrame(results)

# Order columns similar to tu tabla anterior
cols_order = [
    "surface", "method",
    "E0", "E0_L", "E0_U",
    "Emax", "Emax_L", "Emax_U",
    "ED50", "ED50_L", "ED50_U",
    "ED90", "ED90_L", "ED90_U",
    "monotonicity_r", "monotonicity_flag",
]
res_df = res_df.reindex(columns=cols_order)

res_df.to_csv(OUTPUT_TABLE, index=False)
print(f"\nEmax parameter table saved to: {OUTPUT_TABLE}")

# -------------------------------------------------------------------
# 5. Plotting (journal-ready 2x2 layout)
# -------------------------------------------------------------------

# Build mapping for colors/markers if you quieres diferenciarlos
surface_method_pairs = list(zip(df[COL_SURFACE], df[COL_METHOD]))
unique_pairs = sorted(set(surface_method_pairs))

# Pre-define panel layout (like your current Fig. 9)
fig, axes = plt.subplots(2, 2, figsize=(8, 6))
axes = axes.ravel()

# Panel 0: all together
ax_all = axes[0]
ax_all.set_title("Dose–response by surface–method")

for (surface, method) in unique_pairs:
    sub = df[(df[COL_SURFACE] == surface) & (df[COL_METHOD] == method)]
    dose = sub[COL_DOSE].values
    eff = sub[COL_EFF].values

    label = f"{surface} – {method}"

    # scatter points
    ax_all.scatter(dose, eff, label=label, alpha=0.8)

    # overlay Emax curve if we have parameters
    row_match = res_df[(res_df["surface"] == surface) &
                       (res_df["method"] == method)]
    if not row_match.empty:
        d_grid = np.linspace(dose.min(), dose.max(), 200)
        E0_hat = row_match["E0"].iloc[0]
        Emax_hat = row_match["Emax"].iloc[0]
        ED50_hat = row_match["ED50"].iloc[0]
        y_fit = emax_model(d_grid, E0_hat, Emax_hat, ED50_hat)
        ax_all.plot(d_grid, y_fit, linestyle="--")

ax_all.set_xlabel("Dose (mL/L or ppm)")
ax_all.set_ylabel("Efficacy (%)")
ax_all.set_ylim(0, 105)
ax_all.legend(fontsize=7, loc="lower right")

# Panels individuales por combinación (para emular B, C, D)
for idx, (surface, method) in enumerate(unique_pairs, start=1):
    if idx >= len(axes):
        break
    ax = axes[idx]
    sub = df[(df[COL_SURFACE] == surface) & (df[COL_METHOD] == method)]
    dose = sub[COL_DOSE].values
    eff = sub[COL_EFF].values

    ax.scatter(dose, eff, alpha=0.8)
    ax.set_title(f"{surface} – {method}")
    ax.set_xlabel("Dose (mL/L or ppm)")
    ax.set_ylabel("Efficacy (%)")
    ax.set_ylim(0, 105)

    row_match = res_df[(res_df["surface"] == surface) &
                       (res_df["method"] == method)]
    if not row_match.empty:
        d_grid = np.linspace(dose.min(), dose.max(), 200)
        E0_hat = row_match["E0"].iloc[0]
        Emax_hat = row_match["Emax"].iloc[0]
        ED50_hat = row_match["ED50"].iloc[0]
        y_fit = emax_model(d_grid, E0_hat, Emax_hat, ED50_hat)
        ax.plot(d_grid, y_fit, linestyle="--")

plt.tight_layout()
plt.savefig(FIGURE_FILE, dpi=300)
print(f"Figure saved to: {FIGURE_FILE}")
plt.close()
