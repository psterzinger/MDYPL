supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results/new")
out_file <- file.path(results_path, "rescaled-mDYPL-estimates.rda")
n_cores <- 10

library("dplyr")
library("ggplot2")
library("ggpp")
library("patchwork")
library("parallel")
library("RcppNumerical")
## library("brglm2")
devtools::load_all("~/Repositories/brglm2")

source(file.path(supp_path, "code/methods/plot-with-insets.R"))
source(file.path(supp_path, "code/methods/compute-pt.R"))
source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))
source(file.path(supp_path, "code/methods/fit-mdypl.R"))

## Estimates for setting 1
estimate_s1 <- function(kappa, gamma) {
    n <- 1000
    base_beta <- c(-3, -3/2, 0, 3/2, 3)
    p <- n * kappa
    if (p %% 5 > 0) {
        stop("p is not a multiple of 5 for kappa", kappa, ", which is required for the current setting.")
    }
    beta0 <- rep(base_beta, p / 5)
    beta0 <- sort(sqrt(n) * gamma * beta0 / sqrt(sum(beta0^2)))
    X <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
    y <- rbinom(n, 1, plogis(drop(X %*% beta0)))
    coefs <- fit_mdypl(X, y, alpha = 1 / (1 + kappa)) |> coef()
    data.frame(estimate = coefs,
               kappa = kappa,
               gamma = gamma,
               truth = beta0,
               parameter = seq_along(beta0),
               method = "mDYPL")
}

## Estimates for setting 2
estimate_s2 <- function(kappa, gamma) {
    n <- 2000
    base_beta <- c(-10, 10, 0, 0, 0, 0, 0, 0)
    p <- n * kappa
    if (p %% 8 > 0) {
        stop("p is not a multiple of 8 for kappa", kappa, ", which is required for the current setting.")
    }
    beta0 <- rep(base_beta, p / 8)
    beta0 <- sort(sqrt(n) * gamma * beta0 / sqrt(sum(beta0^2)))
    X <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
    y <- rbinom(n, 1, plogis(drop(X %*% beta0)))
    coefs <- fit_mdypl(X, y, alpha = 1 / (1 + kappa)) |> coef()
    data.frame(estimate = coefs,
               kappa = kappa,
               gamma = gamma,
               truth = beta0,
               parameter = seq_along(beta0),
               method = "mDYPL")
}


if (file.exists(out_file)) {
    load(out_file)
} else {
    ## Get phase transition curve
    ns <- 200000
    set.seed(123)
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    ga <- c(0.001, 0.01, seq(0, 20, length = 200))
    pt <- compute_pt(gamma_grid = ga, ncores = n_cores, XZU = xzu)

    ## Get mu, b, sigma for kappa, gamma settings
    kg <- data.frame(kappa = c(0.1, 0.1, 0.2, 0.2, 0.4, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 0.8),
                     gamma = c(10, 18, 6, 14, 2, 10, 18, 6, 14, 2, 10, 18))
    ## Add the setting in Rigon and Aliverti (2023, Section 4.3) and another one where the MLE exists
    kg <- rbind(kg, c(0.2, sqrt(0.9)), c(0.04, 4.5))
    mbs <- matrix(NA, ncol = 3, nrow = nrow(kg))
    se_start <- c(0.5, 1, 1)
    for (i in 1:nrow(kg)) {
        kappa <- kg[i, "kappa"]
        ckappa <- kappa
        gamma <- kg[i, "gamma"]
        mbs[i, ] <- solve_se(kappa, gamma, 1/(1 + kappa), start = se_start)
        if (kappa != ckappa)
            se_start <- se_pars
        cat("kappa =", kappa, "gamma =", gamma, "Done.\n")
    }
    colnames(mbs) <- c("mu", "b", "sigma")
    kgmbs <- cbind(kg, mbs)

    set.seed(123)
    n_reps <- 10
    est_s1 <- est_s2 <- NULL
    for (i in 1:nrow(kgmbs)) {
        kappa <- kgmbs[i, "kappa"]
        gamma <- kgmbs[i, "gamma"]
        mu <- kgmbs[i, "mu"]
        seeds <- generate_unique_seeds(n_reps)
        ests_s1 <- mclapply(1:n_reps, function(i) { set.seed(seeds[i]); estimate_s1(kappa, gamma) }, mc.cores = n_cores)
        ests_s2 <- mclapply(1:n_reps, function(i) { set.seed(seeds[i]); estimate_s2(kappa, gamma) }, mc.cores = n_cores)
        ests_s1 <- do.call("rbind", ests_s1)
        ests_s2 <- do.call("rbind", ests_s2)
        ests_s1$mu <- ests_s2$mu <- mu
        est_s1 <- rbind(est_s1, ests_s1)
        est_s2 <- rbind(est_s2, ests_s2)
        cat("kappa =", kappa, "gamma =", gamma, "Done.\n")
    }

    sum_est_s1 <- est_s1 |> group_by(kappa, gamma, parameter, mu, truth) |>
        summarize(estimate = mean(estimate)) |> data.frame()
    sum_est_s2 <- est_s2 |> group_by(kappa, gamma, parameter, mu, truth) |>
        summarize(estimate = mean(estimate)) |> data.frame()

    save(sum_est_s1, sum_est_s2, pt, file = file.path(results_path, "rescaled-mDYPL-estimates.rda"))
}

## Plots
cols <- hcl.colors(3, palette = "Dark 3")[c(2, 3, 1)]

## Setting 1
p_est_s1 <- plot_with_insets(sum_est_s1, pt, "estimate_and_truth", cols) +
    labs(y = expression(gamma), title = "MDYPL estimator")
p_r_est_s1 <- plot_with_insets(sum_est_s1 |> mutate(estimate = estimate / mu), pt, "estimate_and_truth", cols) +
    labs(y = expression(gamma), title = "rescaled MDYPL estimator")

## Setting 2
p_est_s2 <- plot_with_insets(sum_est_s2, pt, "estimate_and_truth", cols) +
    labs(y = expression(gamma), title = "MDYPL estimator")
p_r_est_s2 <- plot_with_insets(sum_est_s2 |> mutate(estimate = estimate / mu), pt, "estimate_and_truth", cols) +
    labs(y = expression(gamma), title = "rescaled MDYPL estimator")

pdf(file.path(figures_path, paste0("mdypl-vs-truth-1.pdf")), width = 9, height = 4)
print(p_est_s1 + p_r_est_s1)
dev.off()

pdf(file.path(figures_path, paste0("mdypl-vs-truth-2.pdf")), width = 9, height = 4)
print(p_est_s2 + p_r_est_s2)
dev.off()



