# megagravel-thesis
MATLAB code for undergraduate thesis on megagravel transport analysis, UCD 2026
**"Understanding the Movements of Mega Gravels Across the World"**  
David Lee · University College Dublin · 2026  
Supervisor: Prof. Vikram Pakrashi

## Contents

| File | Purpose |
|------|---------|
| `multivariable_regression.m` | Multivariable OLS, log-transformed, and principal component regression. Produces Tables 4-1, 4-2, 4-5 and Figures 4-4 to 4-7. |
| `quantile_regression_tool.m` | Quantile regression pipeline with bootstrap confidence intervals and holdout validation. Produces Tables 4-6 to 4-12 and Figures 4-8 to 4-22. |
| `qreg_lp.m` | Quantile regression solver (linear programming). Called by the pipeline script. |

## Requirements

- MATLAB R2023b or later
- Optimization Toolbox (for `linprog`, used in quantile regression)
- Statistics and Machine Learning Toolbox (for `prctile`, `randsample`)

## Data

The 2013–2014 boulder dataset is not included in this repository. It originates from Cox et al. (2018) and should be obtained from that publication's supplementary material or directly from the authors.

The scripts expect a file named `2014_Boulder_Movements_dataset_of_record.xlsx` in the working directory.

## Running the analysis

1. **Place the dataset** in the same folder as the scripts.

2. **Run `multivariable_regression.m`** to produce the baseline regression results (Section 4.4 of the thesis).

3. **Run `quantile_regression_tool.m` twice** to produce the quantile regression results (Section 4.5):
   - Once with `responseName = 'TotalHorizontalMovement_m__absoluteValue_'` for horizontal movement results
   - Once with `responseName = 'MovementVertical_m_'` for vertical movement results

Set `saveResults = true` at the top of each script to write CSV outputs to disk.

## AI assistance disclosure

AI tools (Claude, Anthropic) were used to assist with code structure, drafting of docstrings, and cleanup of comment style. All algorithmic choices, parameter settings, model specifications, and interpretation of results are my own.

## Licence

MIT — see `LICENSE` file.
