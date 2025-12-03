# Supplementary material for “Diaconis-Ylvisaker prior penalized likelihood for *p*/*n* → *κ* ∈ (0, 1) logistic regression”
Philipp Sterzinger, Ioannis Kosmidis
December 3, 2025

# Directory structure

The directory `code/` contains the R scripts that reproduce all the
numerical results and figures in the manuscript

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

All results are reproducible using R version 4.5.2 (2025-10-31) and the
contributed packages

<table style="width:46%;">
<colgroup>
<col style="width: 26%" />
<col style="width: 19%" />
</colgroup>
<thead>
<tr>
<th>Package</th>
<th>Version</th>
</tr>
</thead>
<tbody>
<tr>
<td>brglm2</td>
<td>1.0.1</td>
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
<td>4.0.1</td>
</tr>
<tr>
<td>ggpp</td>
<td>0.5.9</td>
</tr>
<tr>
<td>glmnet</td>
<td>4.1-10</td>
</tr>
<tr>
<td>memisc</td>
<td>0.99.31.8.3</td>
</tr>
<tr>
<td>mvtnorm</td>
<td>1.3-3</td>
</tr>
<tr>
<td>nleqslv</td>
<td>3.3.5</td>
</tr>
<tr>
<td>patchwork</td>
<td>1.3.2</td>
</tr>
<tr>
<td>RcppNumerical</td>
<td>0.6-0</td>
</tr>
<tr>
<td>tictoc</td>
<td>1.2.1</td>
</tr>
</tbody>
</table>

# Reproducing the results

## Path

All scripts specify the path to the supplementary material path as
`supp_path`. This is currently set to `.` assuming that the working
directory in R is set to the current `git` repository. If this is not
the case for your setup, you should set `supp_path` appropriately.

## Parallel computation

In each of the R scripts relying on parallel computation, `n_cores`
(currently set to `10`) sets the number of cores to use.

Computation in several of the R scripts relies on parallel computing,
which is implemented through the `parallel` R package. The script will
not work on Windows unless `n_cores <- 1` (which will lead in long
compute times and is not recommended) or it is modified to use a
different parallel back-end. All results should be reproducible in
Unix-based systems (e.g. macOS and Linux).

## Details

The following table lists the R scripts that need to be executed in
order to reproduce the results. The table also lists the outputs from
each script, and their label if they are shown in the main text or the
Supplementary Material document. Some of the outputs are intermediate
results, so the scripts should be executed in the order shown.

<table>
<colgroup>
<col style="width: 43%" />
<col style="width: 52%" />
<col style="width: 4%" />
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
href="code/01-rescaled-mDYPL-estimators.R">01-rescaled-mDYPL-estimators.R</a></td>
<td><a
href="rescaled-mDYPL-estimates.rda">rescaled-mDYPL-estimates.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a href="figures/mdypl-vs-truth-1.pdf">mdypl-vs-truth-1.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td></td>
<td><a href="figures/mdypl-vs-truth-2.pdf">mdypl-vs-truth-2.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a
href="code/02-rescaled-PLR-statistics.R">02-rescaled-PLR-statistics.R</a></td>
<td><a
href="results/qqplots-rescaled-plr.rda">qqplots-rescaled-plr.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/qqplots-rescaled-plr.pdf">qqplots-rescaled-plr.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a href="code/03-abias-amse.R">03-abias-amse.R</a></td>
<td><a href="results/abias-amse.rda">abias-amse.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a href="figures/abias-amse.pdf">abias-amse.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a
href="code/04-adaptive-shrinkage-table.R">04-adaptive-shrinkage-table.R</a></td>
<td></td>
<td>Table</td>
</tr>
<tr>
<td><a href="code/05-mu=1.R">05-2-mu=1.R</a></td>
<td><a href="results/mu=1.rda">mu=1.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/alpha-unbiased-min-mse.pdf">alpha-unbiased-min-mse.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a href="code/06-cLS-comparison.R">06-cLS-comparison.R</a></td>
<td><a href="results/">cLS-MDYPL-comparison-*.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/cLS-vs-mDYPL-estimation-inference.pdf">cLS-vs-mDYPL-estimation-inference.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a
href="code/07-min-mse-ml-mdypl-ridge.R">07-min-mse-ml-mdypl-ridge.R</a></td>
<td><a href="results/min_mse.rda">min_mse.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/rlr_mse_scaled_plots.pdf">rlr_mse_scaled_plots.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a
href="code/08-intercept-simul-quantile-tables.R">08-intercept-simul-quantile-tables.R</a></td>
<td><a
href="results/intercept-simu-statistics.rda">intercept-simu-statistics.rda</a></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Table</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Table</td>
</tr>
<tr>
<td><a href="code/09-mfeat-case-study.R">09-mfeat-case-study.R</a></td>
<td><a href="figures/case-study-probs.pdf">case-study-probs.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a href="code/09-mfeat-simulation.R">09-mfeat-simulation.R</a></td>
<td><a href="results/mfeat-simu.rda">mfeat-simu.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/case-study-qqplots.pdf">case-study-qqplots.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a href="code/10-eta-bounds.R">10-eta-bounds.R</a></td>
<td><a href="results/upsilon-gamma.rda">upsilon-gamma.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/cLS-MDYPL-upsilon.pdf">cLS-MDYPL-upsilon.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td></td>
<td><a href="figures/upsilon-gamma.pdf">upsilon-gamma.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td><a
href="11-Bradley-Terry-with-missingness.R">11-Bradley-Terry-with-missingness.R</a></td>
<td><a
href="results/BT-MAR-kappa0.2-gamma-4-alpha-0.7-R-200.rda">BT-MAR-kappa0.2-gamma-4-alpha-0.7-R-200.rda</a></td>
<td></td>
</tr>
<tr>
<td></td>
<td><a
href="figures/bt-centred-estimates.pdf">bt-centred-estimates.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td></td>
<td><a href="figures/bt-inf.pdf">bt-inf.pdf</a></td>
<td>Figure</td>
</tr>
<tr>
<td></td>
<td><a href="figures/bt-qqplots.pdf">bt-qqplots.pdf</a></td>
<td>Figure</td>
</tr>
</tbody>
</table>
