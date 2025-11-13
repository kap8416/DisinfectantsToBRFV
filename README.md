# ToBRFV Disinfectant Efficacy Analysis

This repository contains all scripts and data for analyzing the efficacy of various disinfectants against Tomato brown rugose fruit virus (ToBRFV). The analysis includes lesion-count assays, IRR modeling, and nonlinear Emax dose–response fitting.

## 🧪 Study Overview

The aim of this study is to evaluate the effectiveness of chemical disinfectants in reducing ToBRFV transmission using digital image analysis and statistical modeling. The pipeline integrates image processing, lesion quantification, descriptive statistics, and inferential modeling to identify optimal surface disinfection protocols.

### Experimental Workflow

![Experimental Workflow] (<img width="6000" height="4200" alt="Figure 2-9" src="https://github.com/user-attachments/assets/54a98269-2238-4e82-b514-f839e85e99ec" />)

*Figure: Diagram outlining the experimental workflow used to quantify lesion severity and count from leaf images, perform statistical modeling, and validate treatments using tomato bioassays.*

## 📁 Contents
This repository contains all scripts and data for analyzing the efficacy of various disinfectants against Tomato brown rugose fruit virus (ToBRFV). The analysis includes lesion-count assays, IRR modeling, and nonlinear Emax dose–response fitting.

## 📁 Contents

- `DatosAgronomy.xlsx` — Raw data collected from greenhouse assays.
- `ToBRFV_analysis.py` — Generalized linear modeling (GLM) script to compute IRRs.
- `emax_analysis.py` — Script to clean data, fit dose–response curves (Emax model), and generate plots by surface–method and by product.
- `dose_response_cleaned.csv` — Intermediate CSV generated from raw data.
- `forest_IRR_plastic_spray.png` — Forest plot visualizing IRRs for key conditions.
- `emax_surface_method.png` — Grid plot showing Emax fits across surfaces and methods.
- `emax_by_product_grid.png` — Grid plot showing Emax fits by disinfectant product.

## ▶️ How to Run

### 1. Environment Setup

Ensure Python 3.8+ is installed. Then install required packages:

```bash
pip install pandas numpy matplotlib seaborn scipy openpyxl
2. Run IRR GLM Analysis
bash
Copiar código
python ToBRFV_analysis.py
This script computes Incidence Rate Ratios (IRRs) using negative binomial GLM models and outputs a forest plot.

3. Run Emax Analysis
bash
Copiar código
python emax_analysis.py --input DatosAgronomy.xlsx
This performs data extraction and fits nonlinear dose–response curves (Emax model), generating PNG outputs.

📊 Outputs
emax_surface_method.png — Figure 9A–D. Surface × Method comparison.

emax_by_product_grid.png — Emax fits per product across all surfaces/methods.

forest_IRR_plastic_spray.png — IRR estimates with 95% CI.

📌 Notes
Surface types include: plastic, pruning shears, and hands.

Application methods include: spray and dip.

Doses interpreted in mL/L, ppm, or % depending on label.
