supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "intercept-simu-statistics.rda")
n_cores <- 10

library("brglm2")
library("dplyr")
library("tidyr")
library("RcppNumerical")
library("tictoc")
library("progressr")
handlers("cli")
library("future.apply")
plan(multisession, workers = n_cores)

source(file.path(supp_path, "code/methods/fit-mdypl.R"))

scaled_beta <- function(beta, n, kappa, gamma) {
    p <- floor(n * kappa)
    gamma * beta / sqrt(sum(beta^2) / n)
}

simu_X <- function(n, kappa) {
    p <- floor(n * kappa)
    matrix(rnorm(n * p, 0, 1 / sqrt(n)), n, p)
}

simu_y <- function(X, beta) {
    eta <- drop(X %*% beta)
    rbinom(length(eta), 1, plogis(eta))
}

long2wide <- function(tab) {
    res <- rbind(
        pivot_wider(tab, names_from = "quantile", id_cols = c("kappa", "theta0"),
                    values_from = "z_z"),
        pivot_wider(tab, names_from = "quantile", id_cols = c("kappa", "theta0"),
                    values_from = "z_nz"),
        pivot_wider(tab, names_from = "quantile", id_cols = c("kappa", "theta0"),
                    values_from = "plr_z"))
    res$statistic <- rep(c("z_z", "z_nz", "plr_z"), each = nrow(res) / 3)
    res
}

if (file.exists(out_file)) {
    load(out_file)
} else {
    R <- 5000
    gamma <- sqrt(5)
    n <- 2000
    kt <- expand.grid(theta0 = seq(0.5, 2.5, by = 1),
                      kappa = seq(0.2, 0.8, by = 0.2)) |>
        mutate(p = floor(n * kappa),
               alpha = 1 / (1 + kappa))
    n_settings <- nrow(kt)
    quant <- c(1, 5, 10, 25, 50, 75, 90, 95, 99)
    statistics <- NULL
    for (i in 1:n_settings) {
        theta0 <- kt[i, "theta0"]
        kappa <- kt[i, "kappa"]
        alpha <- kt[i, "alpha"]
        cat(i, "/", n_settings, ": Computing statistics for kappa =", kappa, "and theta0 =", theta0, "\n")
        p <- kt[i, "p"]
        tau = 1 / sqrt(n)
        beta0 <- rep(0, p)
        beta0[sample(1:p, p/4, replace = FALSE)] <- 1
        beta0 <- scaled_beta(beta0, n, kappa, gamma)
        inds_z <- which(beta0 == 0)[1:10]
        coord_z <- inds_z[1]
        coord_nz <- which(beta0 != 0)[1]
        se_consts <- solve_se(kappa, gamma, alpha, intercept = theta0, corrupted = FALSE, start = c(0.5, 1, 1, theta0))
        set.seed(123)
        with_progress(interval = 1, {
            pro <- progressor(R)
            results <- future_lapply(1:R, function(r) {
                X_full <- simu_X(n, kappa)
                X_nest <- X_full[, -inds_z]
                X_full <- cbind(1, X_full)
                X_nest <- cbind(1, X_nest)
                y <- simu_y(X_full, c(theta0, beta0))
                start_full <- c(se_consts[4], beta0 * se_consts[1])
                start_nest <- start_full[-c(inds_z + 1)]
                m_full <- fit_mdypl(X_full, y, alpha, start_full)
                m_nest <- fit_mdypl(X_nest, y, alpha, start_nest)
                plr_z <- 2 * (m_full$pl - m_nest$pl) * se_consts[2] / (kappa * se_consts[3]^2)
                z_z <- sqrt(n) * tau * coef(m_full)[coord_z + 1] / se_consts[3]
                z_nz <- sqrt(n) * tau * (coef(m_full)[coord_nz + 1] - se_consts[1] * beta0[coord_nz]) / se_consts[3]
                pro()
                data.frame(z_z = z_z, z_nz = z_nz, plr_z = plr_z)
            }, future.seed = TRUE)
        })
        results <- do.call("rbind", results)
        tab <- results |>
            reframe(z_z = sapply(qnorm(quant / 100), function(q) 100 * mean(z_z < q)),
                    z_nz = sapply(qnorm(quant / 100), function(q) 100 * mean(z_nz < q)),
                    plr_z = sapply(qchisq(quant / 100, length(inds_z)), function(q) 100 * mean(plr_z < q))) |>
            mutate(quantile = quant, kappa = kappa, theta0 = theta0)
        statistics <- rbind(statistics, tab)
        print(long2wide(tab))
    }
    save(statistics, file = out_file)
}

long2wide(statistics) |>
    filter(statistic == "plr_z")  |>
    data.frame() |> memisc:::toLatex.data.frame(digits = 1)

long2wide(statistics) |>
    filter(statistic == "z_z")  |>
    data.frame() |> memisc:::toLatex.data.frame(digits = 1)

long2wide(statistics) |>
    filter(statistic == "z_nz")  |>
    data.frame() |> memisc:::toLatex.data.frame(digits = 1)
