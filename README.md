# Supplementary material for “Diaconis-Ylvisaker prior penalized likelihood for *p*/*n* → *κ* ∈ (0, 1) logistic regression”

Philipp Sterzinger, Ioannis Kosmidis
December 16, 2024

# Directory structure

The directory `code/` contains the R and Julia scripts that reproduce
all the numerical results and figures in the manuscript

> Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior penalized
> likelihood for *p*/*n* → *κ* ∈ (0, 1) logistic regression.
> https://arxiv.org/abs/2311.07419

and the Supplementary Material document
[`mDYPL-supplementary.pdf`](mDYPL-supplementary.pdf).

The directory `code/methods/` contains methods that the scripts in
`code/` use.

The directory `results/` is populated by files that store the numerical
results the scripts produce.

The directory `figures/` is populated by graphics that the scripts
produce.

# R version and contributed packages

All results are reproducible using R version 4.4.2 (2024-10-31) and the
contributed packages

<table style="width:40%;">
<colgroup>
<col style="width: 26%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr>
<th>Package</th>
<th>Version</th>
</tr>
</thead>
<tbody>
<tr>
<td>colorspace</td>
<td>2.1-1</td>
</tr>
<tr>
<td>detectseparation</td>
<td>0.3</td>
</tr>
<tr>
<td>dplyr</td>
<td>1.1.4</td>
</tr>
<tr>
<td>ggplot2</td>
<td>3.5.1</td>
</tr>
<tr>
<td>ggpp</td>
<td>0.5.8-1</td>
</tr>
<tr>
<td>glmnet</td>
<td>4.1-8</td>
</tr>
<tr>
<td>lmtest</td>
<td>0.9-40</td>
</tr>
<tr>
<td>patchwork</td>
<td>1.3.0</td>
</tr>
<tr>
<td>tidyr</td>
<td>1.3.1</td>
</tr>
</tbody>
</table>

and Julia version 1.10.6 and the contributed packages

<table style="width:46%;">
<colgroup>
<col style="width: 30%" />
<col style="width: 15%" />
</colgroup>
<thead>
<tr>
<th>Package</th>
<th>Version</th>
</tr>
</thead>
<tbody>
<tr>
<td>ColorSchemes</td>
<td>3.27.1</td>
</tr>
<tr>
<td>Cuba</td>
<td>2.3.0</td>
</tr>
<tr>
<td>DataFrames</td>
<td>1.7.0</td>
</tr>
<tr>
<td>Distributions</td>
<td>0.25.114</td>
</tr>
<tr>
<td>FastGaussQuadrature</td>
<td>1.0.2</td>
</tr>
<tr>
<td>FiniteDiff</td>
<td>2.26.2</td>
</tr>
<tr>
<td>InvertedIndices</td>
<td>1.3.1</td>
</tr>
<tr>
<td>JLD2</td>
<td>0.5.10</td>
</tr>
<tr>
<td>LaTeXStrings</td>
<td>1.4.0</td>
</tr>
<tr>
<td>LineSearches</td>
<td>7.3.0</td>
</tr>
<tr>
<td>NLsolve</td>
<td>4.5.1</td>
</tr>
<tr>
<td>NonlinearSolve</td>
<td>4.2.0</td>
</tr>
<tr>
<td>Optim</td>
<td>1.10.0</td>
</tr>
<tr>
<td>Plots</td>
<td>1.40.9</td>
</tr>
<tr>
<td>PrettyTables</td>
<td>2.4.0</td>
</tr>
<tr>
<td>PrettyTables</td>
<td>2.4.0</td>
</tr>
<tr>
<td>ProgressMeter</td>
<td>1.10.2</td>
</tr>
<tr>
<td>RCall</td>
<td>0.14.6</td>
</tr>
<tr>
<td>Roots</td>
<td>2.2.2</td>
</tr>
<tr>
<td>SciMLBase</td>
<td>2.68.1</td>
</tr>
<tr>
<td>StatsBase</td>
<td>0.34.4</td>
</tr>
</tbody>
</table>

# Reproducing the results

## Path

All scripts specify the path to the supplementary material path as
`supp_path`. This is currently set to `.` assuming that the working
directory for Julia and R is set to the current `git` repository. If
this is not the case for your setup, you should set `supp_path`
appropriately.

## Parallel computation

In each of the R and Julia scripts relying on parallel computation,
`n_cores` (currently set to `10`) sets the number of cores to use.

Computation on the R scripts that require it relies on parallel
computing, which is implemented through the `parallel` R package. The
script will not work on Windows unless `n_cores <- 1` (which will lead
in long compute times and is not recommended) or it is modified to use a
different parallel back-end. All results should be reproducible in
Unix-based systems (e.g. macOS and Linux).

Parallel computing in the Julia scripts uses the Julia packages
`Distributed` and `Threads`. See [Parallel
Computing](https://docs.julialang.org/en/v1/manual/parallel-computing/)
in Julia’s documentation for more details, and for setting the number of
threads.

## Details

The following table lists the R and Julia scripts that need to be
executed in order to reproduce the results. The table also lists the
outputs from each script, and their label if they are shown in the main
text or the Supplementary Material document. Some of the outputs are
intermediate results, so the scripts should be executed in the order
shown.

<table>
<colgroup>
<col style="width: 48%" />
<col style="width: 45%" />
<col style="width: 5%" />
</colgroup>
<thead>
<tr>
<th>Script</th>
<th>Output</th>
<th>Label</th>
</tr>
</thead>
<tbody>
<tr>
<td><a
href="code/01-1-rescaled-mDYPL-estimators-1.jl">01-1-rescaled-mDYPL-estimators-1.jl</a></td>
<td><a
href="results/rescaled-mDYPL-estimates-figure-1.rda">rescaled-mDYPL-estimates-figure-1.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/01-1-rescaled-mDYPL-estimators-2.jl">01-1-rescaled-mDYPL-estimators-2.jl</a></td>
<td><a
href="results/rescaled-mDYPL-estimates-figure-2.rda">rescaled-mDYPL-estimates-figure-2.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/01-2-rescaled-mDYPL-estimators-outputs-1.R">01-2-rescaled-mDYPL-estimators-outputs-1.R</a></td>
<td><a href="figures/mdypl-vs-truth-1.pdf">mdypl-vs-truth-1.pdf</a></td>
<td>Figure 1</td>
</tr>
<tr>
<td><a
href="code/01-2-rescaled-mDYPL-estimators-outputs-2.R">01-2-rescaled-mDYPL-estimators-outputs-2.R</a></td>
<td><a href="figures/mdypl-vs-truth-2.pdf">mdypl-vs-truth-2.pdf</a></td>
<td>Figure S1</td>
</tr>
<tr>
<td><a
href="code/02-1-rescaled-PLR-statistics.jl">02-1-rescaled-PLR-statistics.jl</a></td>
<td><a
href="results/qqplots-rescaled-plr.rda">qqplots-rescaled-plr.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/02-2-rescaled-PLR-statistics-outputs.R">02-2-rescaled-PLR-statistics-outputs.R</a></td>
<td><a
href="figures/qqplots-rescaled-plr.pdf">qqplots-rescaled-plr.pdf</a></td>
<td>Figure 2</td>
</tr>
<tr>
<td><a href="code/03-1-abias-amse.jl">03-1-abias-amse.jl</a></td>
<td><a href="results/abias-amse.rda">abias-amse.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/03-2-abias-amse-outputs.R">03-2-abias-amse-outputs.R</a></td>
<td><a href="figures/abias-amse.pdf">abias-amse.pdf</a></td>
<td>Figure 3</td>
</tr>
<tr>
<td><a
href="code/04-adaptive-shrinkage-table.jl">04-adaptive-shrinkage-table.jl</a></td>
<td></td>
<td>Table 1</td>
</tr>
<tr>
<td><a href="code/05-1-min-mse.jl">05-1-min-mse.jl</a></td>
<td><a href="results/alpha-min-mse.rda">alpha-min-mse.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/05-2-mu=1.jl">05-2-mu=1.jl</a></td>
<td><a href="results/mu=1.rda">mu=1.rda</a></td>
<td></td>
</tr>
<tr>
<td><a href="code/05-3-mu=1-outputs.R">05-3-mu=1-outputs.R</a></td>
<td><a
href="figures/alpha-unbiased-min-mse.pdf">alpha-unbiased-min-mse.pdf</a></td>
<td>Figure 4</td>
</tr>
<tr>
<td><a
href="code/06-1-cLS-comparison.jl">06-1-cLS-comparison.jl</a></td>
<td><a href="results/">cLS-MDYPL-comparison-*.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/06-2-cLS-comparison-outputs.R">06-2-cLS-comparison-outputs.R</a></td>
<td><a
href="figures/cLS-vs-mDYPL-estimation-inference.pdf">cLS-vs-mDYPL-estimation-inference.pdf</a></td>
<td>Figure 5</td>
</tr>
<tr>
<td><a
href="code/07-1-min-mse-ml-mdypl-ridge.jl">07-1-min-mse-ml-mdypl-ridge.jl</a></td>
<td><a
href="results/state_evolution_mDYPL.rda">state_evolution_mDYPL.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="results/state_evolution_ridge.rda">state_evolution_ridge.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="results/state_evolution_ML.rda">state_evolution_ML.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/07-2-min-mse-ml-mdypl-ridge-outputs.R">07-2-min-mse-ml-mdypl-ridge-outputs.R</a></td>
<td><a
href="figures/rlr_mse_scaled_plots.pdf">rlr_mse_scaled_plots.pdf</a></td>
<td>Figure 6</td>
</tr>
<tr>
<td><a
href="code/08-1-intercept-simul-quantile-tables.jl">08-1-intercept-simul-quantile-tables.jl</a></td>
<td><a
href="results/intercept_simul-plr-stat-quants.csv">intercept_simul-plr-stat-quants.csv</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="results/intercept-simul-z-stat-quants.csv">intercept-simul-z-stat-quants.csv</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td></td>
<td>Table 2</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Table S1</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Table S2</td>
</tr>
<tr>
<td><a
href="code/09-1-mfeat-case-study.R">09-1-mfeat-case-study.R</a></td>
<td><a href="figures/case-study-probs.pdf">case-study-probs.pdf</a></td>
<td>Figure 7</td>
</tr>
<tr>
<td><a
href="code/09-2-mfeat-simulation.R">09-2-mfeat-simulation.R</a></td>
<td><a href="results/">fou+kar-simu*.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/09-3-mfeat-simulation-outputs.R">09-3-mfeat-simulation-outputs.R</a></td>
<td><a
href="figures/case-study-qqplots.pdf">case-study-qqplots.pdf</a></td>
<td>Figure 7</td>
</tr>
<tr>
<td><a href="code/10-beta-norm-simul.jl">10-beta-norm-simul.jl</a></td>
<td><a href="results/">beta_l2_length*.jld2</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a href="figures/beta_norm_alphas.pdf">beta_norm_alphas.pdf</a></td>
<td>Figure S2</td>
</tr>
<tr>
<td></td>
<td><a href="figures/beta_norm.pdf">beta_norm.pdf</a></td>
<td>Figure S2</td>
</tr>
<tr>
<td><a href="code/11-1-eta-bounds.jl">11-1-eta-bounds.jl</a></td>
<td><a href="results/upsilon-gamma.rda">upsilon-gamma.rda</a></td>
<td></td>
</tr>
<tr>
<td><a
href="code/11-2-eta-gamma-outputs.R">11-2-eta-gamma-outputs.R</a></td>
<td><a
href="figures/cLS-MDYPL-upsilon.pdf">cLS-MDYPL-upsilon.pdf</a></td>
<td>Figure S3</td>
</tr>
<tr>
<td></td>
<td><a href="figures/upsilon-gamma.pdf">upsilon-gamma.pdf</a></td>
<td>Figure 8</td>
</tr>
</tbody>
</table>
