supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
n_cores <- 10

library("dplyr")
library("ggplot2")
library("patchwork")
library("RcppNumerical")
library("brglm2")
library("future.apply")
plan(multisession, workers = n_cores)

source(file.path(supp_path, "code/methods/compute-pt.R"))
source(file.path(supp_path, "code/methods/fit-mdypl.R"))

## Estimates for setting 1
estimate_s <- function(setting) {
    n <- 1000
    if (setting == "a") {
        base_beta <- c(-3, -3/2, 0, 3/2, 3)
        kappa <- 0.2
        gamma <- sqrt(0.9)
    } else if (setting == "b") {
        base_beta <- c(-3, -3/2, 0, 3/2, 3)
        kappa <- 0.05
        gamma <- 10
    } else if (setting == "c") {
        base_beta <- c(-3, -3/2, 0, 3/2, 8)
        kappa <- 0.2
        gamma <- sqrt(3.1)
    }
    p <- n * kappa
    beta0 <- rep(base_beta, p / 5)
    beta0 <- sort(sqrt(n) * gamma * beta0 / sqrt(sum(beta0^2)))
    X <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
    y <- rbinom(n, 1, plogis(drop(X %*% beta0)))
    coefs <- fit_mdypl(X, y, alpha = 1 / (1 + kappa), start = beta0) |> coef()
    data.frame(estimate = coefs,
               kappa = kappa,
               gamma = gamma,
               truth = beta0,
               parameter = seq_along(beta0),
               method = "mDYPL")
}

aggegates <- function(results) {
    kappa <- unique(results$kappa)
    gammasq <- unique(results$gamma)^2
    biases <- with(results, c(tapply(estimate - truth, truth, mean), mean(estimate - truth)))
    mses <- with(results, c(tapply((estimate - truth)^2, truth, mean), mean((estimate - truth)^2)))
    out <- rbind(c(kappa, gammasq, biases),
                 c(kappa, gammasq, sqrt(mses)))
    colnames(out)[c(1:2, 8)] <- c("kappa", "gammasq", "Aggregate")
    colnames(out)[3:7] <- round(as.numeric(colnames(out)[3:7]), 3)
    rownames(out) <- c("bias", "rootMSE")
    out
}

n_reps <- 5000
set.seed(123)
ests_sa <- future_replicate(n_reps, estimate_s("a"), future.seed = TRUE, simplify = FALSE)
ests_sb <- future_replicate(n_reps, estimate_s("b"), future.seed = TRUE, simplify = FALSE)
ests_sc <- future_replicate(n_reps, estimate_s("c"), future.seed = TRUE, simplify = FALSE)
ests_sa <- do.call("rbind", ests_sa)
ests_sb <- do.call("rbind", ests_sb)
ests_sc <- do.call("rbind", ests_sc)

round(aggegates(ests_sa), 3)
round(aggegates(ests_sb), 3)
round(aggegates(ests_sc), 3)
