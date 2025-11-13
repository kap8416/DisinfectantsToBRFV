#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ToBRFV_analysis.py

Full analysis pipeline for ToBRFV disinfectant bioassays:
1. Load and parse Excel file with product blocks
2. Reshape into long format (Product, Treatment, SurfMeth, Dose, Replicate, Lesions)
3. Fit Negative Binomial GLMs per product to obtain IRRs:
   a) Baseline = plastic_spray (if present; otherwise first category)
   b) Baseline = positive_control (POSITIVE WATER)
4. Save IRR tables and generate a forest plot for plastic_spray baseline
"""

import os
import re
import argparse

import numpy as np
import pandas as pd
import statsmodels.api as sm
from patsy import dmatrices
import matplotlib.pyplot as plt


# ---------------------------------------------------------
# 1. Helpers to parse the Excel file
# ---------------------------------------------------------

def load_and_parse_excel(path: str) -> pd.DataFrame:
    """
    Load DatosAgronomy.xlsx and parse it into a long-format DataFrame.

    Expected structure (single sheet):
        - Rows where 'Treatment' is a product name (e.g. 'SOAP', 'GLUTARALDEHÍDO', etc.)
          and lesion columns are NaN.
        - Rows below are actual treatments with lesion counts until the next product
          header or end of file.
        - Columns:
            Treatment, Lesions plant 1, Lesions plant 2, Lesions plant 3, (Average...)

    Returns
    -------
    df_long : DataFrame with columns
        Product, Treatment, Replicate, Lesions
    """
    print(f"Loading Excel file: {path}")
    df_raw = pd.read_excel(path)
    df_raw.columns = [str(c).strip() for c in df_raw.columns]

    # Ensure lesion columns exist with these names
    lesion_cols = [c for c in df_raw.columns if "Lesions" in c]
    if len(lesion_cols) < 3:
        raise ValueError(
            f"Expected at least three lesion columns containing 'Lesions', found: {lesion_cols}"
        )

    # Force 'Treatment' to string
    if "Treatment" not in df_raw.columns:
        raise ValueError("Column 'Treatment' not found in Excel file.")

    # Parse product blocks
    rows_clean = []
    current_product = None

    for _, row in df_raw.iterrows():
        tr = row["Treatment"]
        tr_str = "" if pd.isna(tr) else str(tr).strip()

        # Skip completely empty rows
        if tr_str == "" and all(pd.isna(row[c]) for c in lesion_cols):
            continue

        # Heuristic: if lesion columns are all NaN and Treatment is non-empty,
        # treat this as a product header.
        if all(pd.isna(row[c]) for c in lesion_cols) and tr_str != "":
            current_product = tr_str
            continue

        # Otherwise, we are in a treatment row under current_product
        if current_product is None:
            # If data is messy and starts without a product header, skip
            continue

        # Build record for each replicate
        for r_idx, col in enumerate(lesion_cols, start=1):
            value = row[col]
            if pd.isna(value):
                continue
            rows_clean.append(
                {
                    "Product": current_product,
                    "Treatment": tr_str,
                    "Replicate": r_idx,
                    "Lesions": float(value),
                }
            )

    df_long = pd.DataFrame(rows_clean)
    if df_long.empty:
        raise ValueError("Parsed DataFrame is empty. Check Excel structure.")

    print(f"Parsed {len(df_long)} lesion measurements across products.")
    return df_long


# ---------------------------------------------------------
# 2. Derive SurfMeth, Dose, and control flags
# ---------------------------------------------------------

def parse_surfmeth_and_dose(df: pd.DataFrame) -> pd.DataFrame:
    """
    From the Treatment string, derive:
      - SurfMeth: e.g. 'plastic_spray', 'pruning_shears_dip', 'hands_spray',
                  'healthy_control', 'positive_control'
      - Dose (float, if applicable)

    Assumptions:
      - Strings like:
            '5 mL PLASTIC SPRAY'
            '5 mL PRUNING SHEARS SPRAY'  (o 'TIJERAS ASP')
            '150 ppm HANDS ASP'
            'POSITIVA AGUA' / 'POSITIVE WATER'
            'HEALTHY' / 'SANAS'
    """

    def detect_surface(t: str) -> str:
        t_up = t.upper()
        if "PLAST" in t_up:
            return "plastic"
        if "TIJERAS" in t_up or "PRUNING" in t_up or "SHEARS" in t_up:
            return "pruning_shears"
        if "MANO" in t_up or "HAND" in t_up:
            return "hands"
        if "HEALTHY" in t_up or "SANAS" in t_up:
            return "healthy"
        if "POSITIVE" in t_up or "POSITIVA" in t_up:
            return "positive"
        # Fallback generic
        return "unknown"

    def detect_method(t: str) -> str:
        t_up = t.upper()
        if "SPRAY" in t_up or "ASP" in t_up:
            return "spray"
        if "DIP" in t_up or "INM" in t_up:
            return "dip"
        # Controls or unknown
        return "none"

    def detect_dose(t: str) -> float:
        # Look for first numeric chunk; interpret as float
        matches = re.findall(r"(\d+(\.\d+)?)", t)
        if not matches:
            return np.nan
        return float(matches[0][0])

    surf_list = []
    meth_list = []
    dose_list = []
    surfmeth_list = []

    for tr in df["Treatment"]:
        t = str(tr)
        surf = detect_surface(t)
        meth = detect_method(t)
        dose = detect_dose(t)

        if surf in ["healthy", "positive"]:
            surfmeth = f"{surf}_control"
        else:
            surfmeth = f"{surf}_{meth}"

        surf_list.append(surf)
        meth_list.append(meth)
        dose_list.append(dose)
        surfmeth_list.append(surfmeth)

    df = df.copy()
    df["Surface"] = surf_list
    df["Method"] = meth_list
    df["Dose"] = dose_list
    df["SurfMeth"] = surfmeth_list

    return df


# ---------------------------------------------------------
# 3. Negative Binomial GLM per product
# ---------------------------------------------------------

def fit_nb_glm_by_product(df: pd.DataFrame,
                          baseline_mode: str = "plastic_spray") -> pd.DataFrame:
    """
    Fit a Negative Binomial GLM per product and compute IRRs
    for each SurfMeth level relative to a baseline.

    baseline_mode:
        - "plastic_spray" → internal baseline SurfMeth == 'plastic_spray'
                             (if missing, use first category)
        - "positive_control" → baseline SurfMeth == 'positive_control'
                               (if missing, use first category)

    Returns
    -------
    irr_table : DataFrame with columns
        Product, Baseline, SurfMeth, IRR, LCL, UCL, pvalue, N
    """
    results = []

    products = sorted(df["Product"].unique())
    for prod in products:
        df_prod = df[df["Product"] == prod].copy()

        # Ensure Lesions is numeric
        df_prod["Lesions"] = pd.to_numeric(df_prod["Lesions"], errors="coerce")
        df_prod = df_prod.dropna(subset=["Lesions"])

        if df_prod.empty:
            continue

        # Determine baseline category
        cats = sorted(df_prod["SurfMeth"].unique())
        if not cats:
            continue

        if baseline_mode == "plastic_spray":
            desired_baseline = "plastic_spray"
        elif baseline_mode == "positive_control":
            desired_baseline = "positive_control"
        else:
            raise ValueError("baseline_mode must be 'plastic_spray' or 'positive_control'.")

        if desired_baseline in cats:
            baseline = desired_baseline
        else:
            # Fallback: use first category alphabetically
            baseline = cats[0]

        # Reorder categories so that baseline is the reference
        cats_ordered = [baseline] + [c for c in cats if c != baseline]
        df_prod["SurfMeth"] = pd.Categorical(df_prod["SurfMeth"],
                                             categories=cats_ordered,
                                             ordered=False)

        # Build design matrices
        # Intercept + C(SurfMeth) will use first category as baseline
        formula = "Lesions ~ C(SurfMeth)"
        y, X = dmatrices(formula, data=df_prod, return_type="dataframe")

        # Fit Negative Binomial GLM
        nb_family = sm.families.NegativeBinomial()
        model = sm.GLM(y, X, family=nb_family)
        fit_res = model.fit()

        coefs = fit_res.params
        conf = fit_res.conf_int()
        pvals = fit_res.pvalues

        # Baseline row
        n_baseline = df_prod[df_prod["SurfMeth"] == baseline].shape[0]
        results.append(
            {
                "Product": prod,
                "Baseline": baseline,
                "SurfMeth": baseline,
                "IRR": 1.0,
                "LCL": 1.0,
                "UCL": 1.0,
                "pvalue": np.nan,
                "N": n_baseline,
            }
        )

        # Non-baseline categories
        for cat in cats_ordered:
            if cat == baseline:
                continue

            coef_name = f"C(SurfMeth)[T.{cat}]"
            if coef_name not in coefs.index:
                # This category did not appear in the model (e.g. empty)
                continue

            beta = coefs[coef_name]
            lcl = conf.loc[coef_name, 0]
            ucl = conf.loc[coef_name, 1]
            pval = pvals[coef_name]

            irr = float(np.exp(beta))
            lcl_irr = float(np.exp(lcl))
            ucl_irr = float(np.exp(ucl))

            n_cat = df_prod[df_prod["SurfMeth"] == cat].shape[0]

            results.append(
                {
                    "Product": prod,
                    "Baseline": baseline,
                    "SurfMeth": cat,
                    "IRR": irr,
                    "LCL": lcl_irr,
                    "UCL": ucl_irr,
                    "pvalue": pval,
                    "N": n_cat,
                }
            )

    irr_table = pd.DataFrame(results)
    return irr_table


# ---------------------------------------------------------
# 4. Forest plot for IRR (plastic_spray baseline)
# ---------------------------------------------------------

def plot_forest_irr(irr_table: pd.DataFrame, out_png: str, title: str):
    """
    Generate a forest plot of IRRs (log scale) from the IRR table.
    Significant comparisons (p < 0.05) are colored blue; others gray.
    """
    if irr_table.empty:
        print("IRR table is empty. Skipping forest plot.")
        return

    # Exclude baseline rows from plotting (IRR = 1, p = NaN)
    df_plot = irr_table.copy()
    df_plot = df_plot.sort_values(["Product", "SurfMeth"])

    # Create a combined label for y-axis
    df_plot["Label"] = df_plot["Product"] + " | " + df_plot["SurfMeth"]

    # Sort for a nicer visual (reverse order so first is at top)
    df_plot = df_plot.sort_values("Label", ascending=True)
    df_plot = df_plot.reset_index(drop=True)

    y_pos = np.arange(len(df_plot))
    irr = df_plot["IRR"].values
    lcl = df_plot["LCL"].values
    ucl = df_plot["UCL"].values
    pvals = df_plot["pvalue"].values

    # Significance: p < 0.05 and not NaN
    sig = (pvals < 0.05) & ~np.isnan(pvals)

    colors = np.where(sig, "blue", "gray")

    fig, ax = plt.subplots(figsize=(7, max(6, 0.3 * len(df_plot))))

    # Plot CIs
    ax.hlines(y=y_pos, xmin=lcl, xmax=ucl, color=colors, alpha=0.7)
    # Plot points
    ax.scatter(irr, y_pos, color=colors, zorder=3)

    # Vertical reference line IRR = 1
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1)

    ax.set_xscale("log")
    ax.set_xlabel("Incidence Rate Ratio (IRR, log scale)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot["Label"])
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"Forest plot saved to: {out_png}")


# ---------------------------------------------------------
# 5. Main
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ToBRFV disinfectant analysis: NB-GLM IRRs and forest plots."
    )
    parser.add_argument(
        "--input",
        "-i",
        default="DatosAgronomy.xlsx",
        help="Input Excel file (default: DatosAgronomy.xlsx)",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default="ToBRFV_results",
        help="Output directory (default: ToBRFV_results)",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Load and parse
    df_long = load_and_parse_excel(args.input)

    # 2) Derive SurfMeth and Dose
    df_long = parse_surfmeth_and_dose(df_long)

    # 3a) IRR with plastic_spray baseline
    print("Fitting NB-GLM with baseline = plastic_spray (if present)...")
    irr_plastic = fit_nb_glm_by_product(df_long, baseline_mode="plastic_spray")
    out_csv_plastic = os.path.join(args.outdir, "IRR_plastic_spray_baseline.csv")
    irr_plastic.to_csv(out_csv_plastic, index=False)
    print(f"IRR table (plastic_spray baseline) saved to: {out_csv_plastic}")

    # 3b) IRR with positive_control baseline
    print("Fitting NB-GLM with baseline = positive_control (if present)...")
    irr_positive = fit_nb_glm_by_product(df_long, baseline_mode="positive_control")
    out_csv_positive = os.path.join(args.outdir, "IRR_positive_control_baseline.csv")
    irr_positive.to_csv(out_csv_positive, index=False)
    print(f"IRR table (positive_control baseline) saved to: {out_csv_positive}")

    # 4) Forest plot for plastic_spray baseline
    out_png = os.path.join(args.outdir, "forest_IRR_plastic_spray.png")
    plot_forest_irr(
        irr_table=irr_plastic,
        out_png=out_png,
        title="IRR by product and surface–method (baseline: plastic_spray)",
    )

    print("Done.")


if __name__ == "__main__":
    main()
