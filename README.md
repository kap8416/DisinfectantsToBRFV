# ToBRFV Disinfectant Efficacy Analysis

This repository contains all scripts and data for analyzing the efficacy of various disinfectants against Tomato brown rugose fruit virus (ToBRFV). The analysis includes lesion-count assays, IRR modeling, and nonlinear Emax dose–response fitting.

## 🧪 Study Overview

The aim of this study is to evaluate the effectiveness of chemical disinfectants in reducing ToBRFV transmission using digital image analysis and statistical modeling. The pipeline integrates image processing, lesion quantification, descriptive statistics, and inferential modeling to identify optimal surface disinfection protocols.

### Experimental Workflow

<img width="6000" height="4200" alt="Figure 2-9" src="https://github.com/user-attachments/assets/54a98269-2238-4e82-b514-f839e85e99ec" />

*Figure: Diagram outlining the experimental workflow used to quantify lesion severity and count from leaf images, perform statistical modeling, and validate treatments using tomato bioassays.*

## Descriptive Evaluation Prior to Model Fitting
A descriptive statistical analysis was first performed to identify patterns of virucidal performance across disinfectants, doses, surfaces, and application methods. For each treatment, lesion count and leaf damage severity were summarized using means, standard deviations, and 95% confidence intervals. These descriptive metrics enabled the identification of the most and least effective dose–method combinations prior to model-based inference. All calculations and visualizations were conducted in Python and R, using the pandas, numpy, and ggplot2 libraries. This descriptive evaluation provided the empirical foundation for subsequent inferential modeling.

## Statistical Modeling of Surface–Method Interactions

Lesion count data were statistically evaluated to determine the effects of disinfectant formulation, surface type, and application method. Data exploration revealed overdispersion, leading to the use of a Negative Binomial Generalized Linear Model (NB-GLM) as the primary analytical framework. Model selection was based on residual diagnostics and Akaike Information Criterion, with the NB-GLM outperforming Poisson and Quasi-Poisson alternatives. The model included disinfectant, surface, and method as fixed effects, along with their pairwise interactions. Incidence Rate Ratios (IRRs) were calculated relative to the plastic–spray reference group to estimate the relative infection risk. Confidence intervals (95%) and Wald tests (p < 0.05) were used to assess significance. For disinfectants showing a monotonic response to concentration, a nonlinear Emax model was applied to characterize dose-dependent efficacy. Parameter estimation was performed using weighted least squares and Levenberg–Marquardt optimization, with model performance evaluated through residual analysis and AIC.

## 📁 Contents
This repository contains all scripts and data for analyzing the efficacy of various disinfectants against Tomato brown rugose fruit virus (ToBRFV). The analysis includes lesion-count assays, IRR modeling, and nonlinear Emax dose–response fitting.

## 📁 

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
python ToBRFV_analysis.py
This script computes Incidence Rate Ratios (IRRs) using negative binomial GLM models and outputs a forest plot.

3. Run Emax Analysis
bash
python emax_analysis.py --input DatosAgronomy.xlsx
This performs data extraction and fits nonlinear dose–response curves (Emax model), generating PNG outputs.

📊 Outputs
emax_surface_method.png — Surface × Method comparison.

emax_by_product_grid.png — Emax fits per product across all surfaces/methods.

forest_IRR_plastic_spray.png — IRR estimates with 95% CI.

📌 Notes
Surface types include: plastic, pruning shears, and hands.

Application methods include: spray and dip.

Doses interpreted in mL/L, ppm, or % depending on label.
