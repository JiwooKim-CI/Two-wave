# Complete Two-Wave Mediation Analyses: Seven- and Nine-Point Grids

This package runs both confounding models under three scenarios each: unrestricted, common trend, and ignorability. It enumerates all combinations of the free grid coefficients, saves numerical summaries, produces half-eye and bias-versus-violation figures, and compares the seven- and nine-point results.

## Files

- `simulation_model1.R`: shared-confounder model, three scenarios, standalone.
- `simulation_model2.R`: time-varying-confounder model, three scenarios, standalone.
- `run_7point.R`: run both models on the seven-point grid.
- `run_9point.R`: run both models on the nine-point grid.
- `run_all.R`: run both grids and create comparison outputs.
- `run_grid.R`: launch each model in a separate R process.
- `compare_grids.R`: compare completed seven-/nine-point results.
- `numerical_checks.csv`: independent numerical verification results for this version.

Keep these files in one directory and use that directory as the working directory.

## Install packages once

```r
install.packages(c("data.table", "ggplot2", "ggdist", "jtools", "ragg"))
```

## Run all analyses

In RStudio:

```r
source("run_all.R")
```

On an Azure Linux VM:

```sh
Rscript --vanilla run_all.R
```

The runner executes the two seven-point analyses, then the two nine-point analyses, then the comparison. Each analysis runs in its own R process. Results are organized as:

```text
results/
  grid7/model1/
  grid7/model2/
  grid9/model1/
  grid9/model2/
  grid_comparison/
```

To run a single grid:

```r
source("run_7point.R")
# or
source("run_9point.R")
```

After both grids finish, regenerate the comparison with:

```r
source("compare_grids.R")
compare_grids("results")
```

The two model files can also be run independently:

```sh
N_POINTS=9 OUTPUT_DIR=results/grid9/model1 Rscript --vanilla simulation_model1.R
N_POINTS=9 OUTPUT_DIR=results/grid9/model2 Rscript --vanilla simulation_model2.R
```

## Small end-to-end test

```r
source("run_grid.R")
run_grid(3L, root="smoke_results")
```

The three-point grid uses -0.1, 0, and 0.1 and is for testing the workflow. The production runners explicitly select seven or nine points.

## Parameters and enumeration

The seven-point grid contains 0, ±0.1, ±sqrt(0.1), and ±fourth-root(0.1). The nine-point grid adds ±eighth-root(0.1). Beta traverses all grid values in both models. Alpha and gamma are fixed at 0.20.

| Independently generated scenario | Seven-point candidates | Nine-point candidates |
|---|---:|---:|
| Shared unrestricted | 282,475,249 | 3,486,784,401 |
| Time-varying unrestricted | 1,977,326,743 | 31,381,059,609 |
| Shared common trend | 823,543 | 4,782,969 |
| Time-varying common trend | 5,764,801 | 43,046,721 |

Ignorability results are obtained from each unrestricted grid using the model-specific delta3*delta4=0 condition.

This is an evaluation of population quantities, not repeated generation of finite-sample datasets. Common-trend coefficients are derived before admissibility checks and can fall between or outside the original grid points.

## Corrected common-trend scenarios

### Shared confounder

```r
p_m2_y1 = 0
p_y2_u  = (1 - p_y2_y1) * p_y1_u
p_y2_m1 = (1 - p_y2_y1) * p_y1_m1
```

### Time-varying confounders

The second assignment is evaluated before the third:

```r
p_m2_y1 = 0
p_y2_u2 = p_u2_u1 * (1 - p_y2_y1) * p_y1_u1 *
          (1 - p_m1_u1^2) / (1 - p_u2_u1^2 * p_m1_u1^2)
p_y2_m1 = (1 - p_y2_y1) * (p_y1_m1 + p_m1_u1 * p_y1_u1) -
          p_m1_u1 * p_u2_u1 * p_y2_u2
```

These restrictions make the outcome-only estimator Change(D) unbiased under the maintained structural assumptions. At eta=1, they reduce to the shared-confounder restrictions. Each common-trend chunk checks that the maximum absolute Change(D) bias is at most 1e-8; each ignorability chunk checks the corresponding ANCOVA result.

## Calculation corrections

- Duplicated residual-variance cross terms are corrected before taking square roots.
- Required residual and conditional variances must exceed 1e-10.
- Invalid change-score correlation filters are removed; raw difference-score covariance/variance formulas are used.
- All beta values are enumerated in both models.
- Undefined objects, nonexistent columns, positional indexing, and duplicate exports are removed.
- Change(A) is computed as twice Change(F) after Change(F) exists.
- Best-estimator comparisons use the named, unrounded bias columns, with a 1e-10 tie tolerance.
- CSV frequencies include both bias and the magnitude of the ignorability violation.
- Both common-trend scenarios are constructed directly rather than using the old equality filters.

## Memory and execution controls

The original `CJ()` enumeration and parameter names are preserved. Nested loops divide the grid into chunks; plotting frequencies and numerical sums are accumulated instead of retaining billions of raw rows.

| Environment variable | Default | Purpose |
|---|---|---|
| `N_POINTS` | 7 | Standalone model grid size |
| `INNER_DIM` | 6 | Maximum free dimensions in one inner CJ grid |
| `OUTPUT_DIR` | Timestamped model directory | Standalone result directory |
| `DRAW_PLOTS` | TRUE | Generate figures after calculation |
| `R_THREADS` | 1 | data.table threads within an R process |

For example, to reduce chunk size for the nine-point run:

```sh
INNER_DIM=5 R_THREADS=4 Rscript --vanilla run_9point.R
```

At INNER_DIM=6, a chunk has at most 117,649 rows for seven points or 531,441 for nine. INNER_DIM=5 reduces those numbers to 16,807 and 59,049. Accumulated frequency tables also occupy memory. Set resources according to the target VM.

## Checkpoint and resume

After each completed beta block, `checkpoint.rds` saves accumulated results and the loop position. The checkpoint is written to a temporary file and then renamed. To resume, run the same command with the same output path, grid, and INNER_DIM. Completed beta blocks are skipped; an interrupted block is recomputed. The checkpoint checks the model, grid settings, and calculation/accumulation function definitions.

The production runners use stable result paths so that rerunning `run_all.R` resumes previous work. Re-running a finished analysis regenerates final exports and figures from its checkpoint. For a separate fresh analysis, choose a new output root through `run_grid()` or a new standalone OUTPUT_DIR.

## Numerical outputs

Each model/grid folder contains:

- `bias_summary.csv`: N, mean bias, mean absolute bias, maximum absolute bias, and accumulated sums, calculated before rounding.
- `counts.csv`: candidates and retained counts for general/common-trend generation; ignorability N is in the summary.
- `best_estimator.csv`: counts and proportions for the smallest absolute bias, including ties.
- `violation_curve.csv`: unrounded absolute biases averaged within violation levels (the grouping variable is rounded to 12 decimals).
- `plot_results.rds`: plotting frequencies, summaries, and curve data.
- Original-style per-estimator CSVs containing value, diff, n, and perc. Bias and diff are rounded to four decimals for plotting frequencies; perc is a proportion.
- `checkpoint.rds`: resumable accumulated state.
- `session_info.txt`: grid and R environment.

Time-varying common-trend CSVs use the restored `_2_1_c.csv` suffix, for example `ancova_2_1_c.csv`.

## Figures

The top-to-bottom method order is Naive, ANCOVA, Change (F), Change (A), Change (D).

The half-eye figures use the supplied horizontal `stat_halfeye()` plus boxplot layout, method colors, frequency weights, and `theme_apa()`. Common-trend bandwidth adjustment is 5; other scenarios use 0.5. The viewing window is -0.5 to 0.5 without deleting values before density estimation. PDF size is six by six inches; PNGs are generated when ragg is available.

| Model/scenario | Figure basename |
|---|---|
| Shared unrestricted | new_sim1 |
| Shared common trend | new_sim2 |
| Shared ignorability | new_sim3 |
| Time-varying unrestricted | new_sim2_1 |
| Time-varying common trend | new_sim2_2 |
| Time-varying ignorability | new_sim2_3 |
| Shared bias versus violation | sim_1_ig_2 |
| Time-varying bias versus violation | sim_2_ig_2 |

Line figures are generated automatically from `violation_curve.csv` quantities, using exact accumulated mean absolute biases, distinct black line types, and direct labels. Their axes include the computed curves rather than imposing the old 0.11 upper bound.

## Seven- versus nine-point comparison

`run_all.R` generates:

- `grid_comparison/comparison_7_vs_9.csv`: mean absolute biases, ranks, maximum absolute biases, counts, and differences for every estimator in all six scenarios.
- `grid_comparison/identification_checks.csv`: zero-bias checks for Change(D) under common trend and ANCOVA under ignorability, for both grids/models.
- `grid_comparison/comparison_7_vs_9.pdf` and PNG: a visual comparison.

The nine-point grid expands the coefficient range as well as adding grid values. The comparison documents the numerical changes rather than automatically declaring the conclusions unchanged.

## Verification

This revision was checked by independently constructing structural covariance matrices and regression coefficients for 1,000 candidate configurations in each of eight settings: two models, two grids, and general/common-trend scenarios. Residual variances, retained configurations, and all five estimates agreed to floating-point precision. The numerical results are recorded in `numerical_checks.csv`.

Small full-grid runs verified all six scenarios and the eight figure exports. Production seven- and nine-point results are generated by the supplied runners; the bundled numerical-check file is validation evidence, not the completed production analysis.
