#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Improved and modularized script to:
1. Extract and clean efficacy data from Excel
2. Fit Emax models per surface-method and per product
3. Generate and save all plots to output directory
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

# --- Configuration ---
EXCEL_FILE = "DatosAgronomy.xlsx"
OUTPUT_DIR = "outputs"
CLEAN_CSV = os.path.join(OUTPUT_DIR, "dose_response_cleaned.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Emax Model Function ---
def emax_model(dose, Emax, EC50):
    return Emax * dose / (EC50 + dose)

# --- Step 1: Extract and Clean Data ---
def extract_and_clean_data(excel_file):
    df_all = pd.read_excel(excel_file, sheet_name=None)
    records = []

    for sheet_name, sheet_df in df_all.items():
        sheet_df = sheet_df.dropna(how="all")
        sheet_df.columns = sheet_df.columns.str.strip()
        product = sheet_name.strip().lower()

        pos_control = sheet_df[sheet_df.iloc[:, 0].astype(str).str.upper() == "POSITIVE WATER"]
        if pos_control.empty:
            continue

        pos_mean = float(pos_control.iloc[0, -1])

        for _, row in sheet_df.iterrows():
            treat = str(row[0]).strip().upper()
            if "HEALTHY" in treat or "POSITIVE WATER" in treat:
                continue

            lesions = row[-1]
            if pd.isna(lesions) or pos_mean == 0:
                continue

            if "PLASTIC" in treat:
                surface = "plastic"
            elif "PRUNING SHEARS" in treat:
                surface = "pruning shears"
            elif "HAND" in treat:
                surface = "hand"
            else:
                continue

            if "SPRAY" in treat:
                method = "spray"
            elif "DIP" in treat:
                method = "dip"
            else:
                continue

            match = re.search(r"([\d.]+)\s*(PPM|ML|%)", treat)
            if match:
                dose = float(match.group(1))
            else:
                continue

            efficacy = 100 * (1 - float(lesions) / pos_mean)
            efficacy = max(0, min(100, efficacy))

            records.append({
                "Product": product,
                "Surface": surface,
                "Method": method,
                "Dose": dose,
                "Efficacy": efficacy
            })

    df_clean = pd.DataFrame(records)
    df_clean.to_csv(CLEAN_CSV, index=False)
    print(f"✅ Exported: {CLEAN_CSV}")
    return df_clean

# --- Step 2A: Plot Grid by Surface–Method ---
def plot_by_surface_method(df):
    combinations = [
        ("plastic", "spray"),
        ("pruning shears", "dip"),
        ("pruning shears", "spray"),
    ]
    palette = sns.color_palette("colorblind", n_colors=len(combinations))
    color_map = {combo: palette[i] for i, combo in enumerate(combinations)}

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    ax_combined = axes[0, 0]
    ax_map = {
        ("plastic", "spray"): axes[0, 1],
        ("pruning shears", "dip"): axes[1, 0],
        ("pruning shears", "spray"): axes[1, 1],
    }

    for (surface, method) in combinations:
        sub_df = df[(df.Surface == surface) & (df.Method == method)]
        if sub_df.empty:
            continue

        x, y = sub_df["Dose"], sub_df["Efficacy"]
        color = color_map[(surface, method)]
        Emax_g, EC50_g = y.max(), np.median(x[x > 0])

        try:
            popt, _ = curve_fit(emax_model, x, y, p0=[Emax_g, EC50_g], maxfev=10000)
        except RuntimeError:
            popt = [Emax_g, EC50_g]

        x_curve = np.linspace(x.min(), x.max(), 100)
        y_curve = emax_model(x_curve, *popt)

        ax_combined.scatter(x, y, color=color, label=f"{surface} – {method}")
        ax_combined.plot(x_curve, y_curve, color=color, linestyle="--")

        ax = ax_map[(surface, method)]
        ax.scatter(x, y, color=color)
        ax.plot(x_curve, y_curve, color=color, linestyle="--")
        ax.set_title(f"{surface.capitalize()} – {method.capitalize()}")
        ax.set_ylim(0, 105)

    ax_combined.set_title("Dose–response by surface–method")
    ax_combined.set_ylim(0, 105)
    ax_combined.legend(title="Surface / Method")

    for ax in [ax_combined] + list(ax_map.values()):
        ax.set_xlabel("Dose (mL/L or ppm)")
        ax.set_ylabel("Efficacy (%)")

    plt.savefig(os.path.join(OUTPUT_DIR, "emax_surface_method.png"), dpi=300)
    plt.close()

# --- Step 2B: Plot Grid by Product ---
def plot_by_product(df):
    products = sorted(df["Product"].unique())
    n_cols = 3
    n_rows = int(np.ceil(len(products) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.5 * n_cols, 5 * n_rows), constrained_layout=True)
    axes = axes.flatten()
    palette = sns.color_palette("colorblind", 10)

    for i, product in enumerate(products):
        ax = axes[i]
        prod_df = df[df.Product == product]
        combos = prod_df.groupby(["Surface", "Method"])

        for j, ((surface, method), group) in enumerate(combos):
            x, y = group["Dose"], group["Efficacy"]
            color = palette[j % len(palette)]
            label = f"{surface} – {method}"
            Emax_g, EC50_g = y.max(), np.median(x[x > 0])

            try:
                popt, _ = curve_fit(emax_model, x, y, p0=[Emax_g, EC50_g], maxfev=10000)
            except RuntimeError:
                popt = [Emax_g, EC50_g]

            x_curve = np.linspace(x.min(), x.max() * 1.1, 100)
            y_curve = emax_model(x_curve, *popt)

            ax.scatter(x, y, color=color, edgecolor="k")
            ax.plot(x_curve, y_curve, color=color, linestyle="--", label=label)

        ax.set_title(product.upper(), fontsize=11)
        ax.set_xlabel("Dose (mL/L or ppm)")
        ax.set_ylabel("Efficacy (%)")
        ax.set_ylim(0, 105)
        ax.legend(fontsize=8, loc="lower right")

    for k in range(len(products), len(axes)):
        fig.delaxes(axes[k])

    plt.savefig(os.path.join(OUTPUT_DIR, "emax_by_product_grid.png"), dpi=300)
    plt.close()

# --- Main Routine ---
if __name__ == "__main__":
    df_clean = extract_and_clean_data(EXCEL_FILE)
    plot_by_surface_method(df_clean)
    plot_by_product(df_clean)
    print("✅ All figures saved to outputs/")
