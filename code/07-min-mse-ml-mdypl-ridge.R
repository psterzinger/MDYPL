supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "min_mse.rda")
n_cores <- 10

library("statmod")
library("dplyr")
library("tidyr")
library("ggplot2")
library("patchwork")
library("brglm2")
library("progressr")
handlers("cli")
library("future.apply")
plan(multisession, workers = n_cores)


source(file.path(supp_path, "code/methods/compute-pt.R"))
source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))

alpha_min_mse <- function(start, kappa, gamma, int = c(0, 1), init_iter = 0, gh = gauss.quad(200, "hermite")) {
    g <- function(alpha) {
        pars <- solve_se(kappa, gamma, alpha, start = start, init_iter = init_iter,
                         control = list(ftol = 1e-12), gh = gh)
        sqrt(kappa) * pars[3] / pars[1]
    }
    res <- optimise(g, int, tol = 1e-08)
    c(alpha = res$minimum, rmse = res$objective)
}

lambda_min_mse <- function(start, kappa, gamma, int = c(0, 1), init_iter = 0, gh = gauss.quad(200, "hermite")) {
    g <- function(lambda) {
        pars <- brglm2:::solve_se_ridge(kappa, gamma, lambda, start = start, init_iter = init_iter,
                                        control = list(ftol = 1e-12), gh = gh)
        sqrt(kappa) * pars[3] / pars[1]
    }
    res <- optimise(g, int, tol = 1e-08)
    c(lambda = res$minimum, rmse = res$objective)
}

if (file.exists(out_file)) {
    load(out_file)
} else {
    kappas <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    gammas <- c(1.0, 2.5, 5, 7.5, 10, 12.5, 15.0)
    kg <- expand.grid(kappa = kappas, gamma = gammas)
    n_kg <- nrow(kg)
    ## Compute phase transition
    kg <- as.data.frame(kg)
    colnames(kg) <- c("kappa", "gamma")
    set.seed(123)
    ns <- 200000
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    kg$kappa_diff <- future_sapply(seq.int(nrow(kg)), function(j) {
        gamma <- kg$gamma[j]
        kappa <- kg$kappa[j]
        compute_pt(beta0 = 0, gamma, ncores = 1, XZU =  xzu)$kappa - kappa
    })
    mle_exists <- kg$kappa_diff > 0
    grid_size <- 20
    alphas <- c(seq(0.01, 0.9, length.out = 2 * grid_size),
                seq(0.901, 0.99, length.out = 3 * grid_size))
    lambdas <- c(seq(4.5, 0.1, length.out = 2 * grid_size),
                 seq(0.099, 0.001, length.out = 3 * grid_size))
    n_a <- length(alphas)
    n_l <- length(lambdas)
    ## ML
    start0 <- c(0.7, 2, 2)
    with_progress({
        pro <- progressor(n_kg)
        df_ml <- future_lapply(1:n_kg, function(i) {
            kappa = kg[i, "kappa"]
            gamma = kg[i, "gamma"]
            if (mle_exists[i]) {
                consts <- solve_se(kappa, gamma, 1, start = start0)
                linf <- max(attr(consts, "funcs"))
            } else {
                consts <- rep(NA, 3)
                linf <- NA
            }
            pro(message = paste(linf))
            data.frame(kappa = kappa, gamma = gamma, alpha = 1,
                       mu = consts[1], b = consts[2], sigma = consts[3], linf = linf)
        })
    })
    df_ml <- do.call("rbind", df_ml)
    df_ml$method <- "ML"
    df_ml_min <- df_ml |>
        mutate(mse = kappa * sigma^2 / mu^2)

    ## MDYPL
    ## Get solutions at alphas[1]
    start_sol <- matrix(NA, n_kg, 4)
    start_mdypl <- c(0.8, 0.1, 0.1)
    for (i in 1:n_kg) {
        kappa = kg[i, "kappa"]
        gamma = kg[i, "gamma"]
        alpha <- alphas[1]
        if ((i - 1) %% length(kappas) == 0) {
            start <- if (i == 1) start_mdypl else start_sol[i - length(kappas), 1:3]
        } else {
            start <- start_sol[i - 1, 1:3]
        }
        res <- solve_se(kappa, gamma, alpha,
                        start = start,
                        init_iter = 0, control = list(ftol = 1e-12))
        start_sol[i, 1:3] <- res
        start_sol[i, 4] <- max(abs(attr(res, "funcs")))
        cat(i, "/", n_kg, "kappa =", kappa, "gamma =", gamma, "alpha =", alpha,
            "linf", max(start_sol[, 4], na.rm = TRUE), "\n")
    }
    colnames(start_sol) <- c("mu", "b", "sigma", "linf")
    ## Get solutions for all alphas
    with_progress({
        every <- 20
        pro <- progressor(n_kg * n_a / every)
        results <- future_lapply(1:n_kg, function(i) {
            kappa <- kg[i, "kappa"]
            gamma <- kg[i, "gamma"]
            consts <- matrix(NA, n_a, 4)
            colnames(consts) <- c("mu", "b", "sigma", "linf")
            start <- start_sol[i, c("mu", "b", "sigma")]
            for (ia in 1:n_a) {
                alpha <- alphas[ia]
                start <- if (ia > 1) consts[ia - 1, 1:3] else start
                res <- solve_se(kappa, gamma, alpha, start = start, init_iter = 0,
                                control = list(ftol = 1e-12))
                consts[ia, 1:3] <- res
                consts[ia, 4] <- max(abs(attr(res, "funcs")))
                if (ia %% every == 0)
                    pro(paste("kappa =", kappa, "gamma =", gamma, "alpha =", alpha,
                              "linf", max(consts[, 4], na.rm = TRUE)))
            }
            data.frame(kappa = kappa, gamma = gamma, alpha = alphas, consts)
        })
    })
    df_mdypl <- do.call("rbind", results)
    df_mdypl$method <- "MDYPL"
    ## Refine with more quad points
    df_mdypl_min <- df_mdypl |>
        mutate(mse = kappa * sigma^2 / mu^2) |>
        group_by(kappa, gamma) |>
        slice_min(order_by = mse, with_ties = FALSE) |>
        ungroup()
    with_progress({
        pro <- progressor(nrow(df_mdypl_min))
        alpha_exact <-  future_lapply(1:nrow(df_mdypl_min), function(i) {
            cpars <- df_mdypl_min[i, ]
            kappa <- cpars$kappa
            gamma <- cpars$gamma
            alpha <- cpars$alpha
            mu <- cpars$mu
            b <- cpars$b
            sigma <- cpars$sigma
            int <- alpha + c(-0.1, 0.1)
            int[int < 0] <- 0.001
            int[int > 1] <- 0.999
            out <- alpha_min_mse(c(mu, b, sigma), kappa, gamma, int = int,
                                 init_iter = 0, gh = gauss.quad(1000, "hermite"))
            pro(paste(round(out[1], 3), round(out[2]^2, 3)))
            out
        })
    })
    alpha_exact <- do.call("rbind", alpha_exact)
    df_mdypl_min$alpha <- alpha_exact[, "alpha"]
    df_mdypl_min$mse <- alpha_exact[, "rmse"]^2

    ## Logistic ridge
    ## Get solutions at lambda[1]
    start_sol <- matrix(NA, n_kg, 4)
    start_ridge <- c(0.8, 0.1, 0.1)
    for (i in 1:n_kg) {
        kappa = kg[i, "kappa"]
        gamma = kg[i, "gamma"]
        lambda <- lambdas[1]
        if ((i - 1) %% length(kappas) == 0) {
            start <- if (i == 1) start_ridge else start_sol[i - length(kappas), 1:3]
        } else {
            start <- start_sol[i - 1, 1:3]
        }
        res <- brglm2:::solve_se_ridge(kappa, gamma, lambda, start = start_mdypl,
                                       control = list(ftol = 1e-12))
        start_sol[i, 1:3] <- res
        start_sol[i, 4] <- max(abs(attr(res, "funcs")))
        cat(i, "/", n_kg, "kappa =", kappa, "gamma =", gamma, "lambda =", lambda,
            "linf", max(start_sol[, 4], na.rm = TRUE), "\n")
    }
    colnames(start_sol) <- c("mu", "b", "sigma", "linf")
    ## Get solutions for all lambdas
    with_progress({
        every <- 20
        pro <- progressor(n_kg * n_l / every)
        results <- future_lapply(1:n_kg, function(i) {
            kappa <- kg[i, "kappa"]
            gamma <- kg[i, "gamma"]
            consts <- matrix(NA, n_a, 4)
            colnames(consts) <- c("mu", "b", "sigma", "linf")
            start <- start_sol[i, c("mu", "b", "sigma")]
            for (il in 1:n_l) {
                lambda <- lambdas[il]
                start <- if (il == 1) start else consts[il - 1, 1:3]
                res <- brglm2:::solve_se_ridge(kappa, gamma, lambda, start = start, init_iter = 0,
                                               control = list(ftol = 1e-14))
                consts[il, 1:3] <- res
                consts[il, 4] <- max(abs(attr(res, "funcs")))
                if (il %% every == 0)
                    pro(message = paste("kappa =", kappa, "gamma =", gamma, "lambda =", lambda,
                                        "linf", max(consts[, 4], na.rm = TRUE)))
        }
        data.frame(kappa = kappa, gamma = gamma, lambda = lambdas, consts)
        })
    })
    df_ridge <- do.call("rbind", results)
    df_ridge$method <- "ridge"
    ## Refine with more quad points
    df_ridge_min <- df_ridge |>
        mutate(mse = kappa * sigma^2 / mu^2) |>
        group_by(kappa, gamma) |>
        slice_min(order_by = mse, with_ties = FALSE) |>
        ungroup()
    with_progress({
        pro <- progressor(nrow(df_ridge_min))
        lambda_exact <-  future_lapply(1:nrow(df_ridge_min), function(i) {
            cpars <- df_ridge_min[i, ]
            kappa <- cpars$kappa
            gamma <- cpars$gamma
            lambda <- cpars$lambda
            mu <- cpars$mu
            b <- cpars$b
            sigma <- cpars$sigma
            int <- lambda + c(-0.5, 0.5)
            int[int < 0] <- 0.0001
            out <- lambda_min_mse(c(mu, b, sigma), kappa, gamma, int = int,
                                  init_iter = 0, gh = gauss.quad(1000, "hermite"))
            pro(paste(round(out[1], 3), round(out[2]^2, 3)))
            out
        })
    })
    lambda_exact <- do.call("rbind", lambda_exact)
    df_ridge_min$lambda <- lambda_exact[, "lambda"]
    df_ridge_min$mse <- lambda_exact[, "rmse"]^2
    ## Compute MSEs and find minimum MSE
    all_mses <- rbind(df_ml_min[c("kappa", "gamma", "mse", "method")],
                      df_mdypl_min[c("kappa", "gamma", "mse", "method")],
                      df_ridge_min[c("kappa", "gamma", "mse", "method")])
    all_mses <- all_mses |>
        pivot_wider(names_from = method, values_from = mse) |>
        mutate(ratio_ridge = MDYPL / ridge,
               ratio_ML = MDYPL / ML,
               mle_exists = !is.na(ratio_ML)) |>
        data.frame()
    all_mses <- all_mses |> mutate(gamma_lab = factor(paste("gamma ==", gamma),
                                                      levels = paste("gamma ==", unique(gamma)),
                                                      ordered = TRUE))
    save(df_ml_min, df_mdypl_min, df_ridge_min, all_mses, file = out_file)
}

mDYPL_ridge <- ggplot(all_mses) +
    geom_col(aes(kappa, sqrt(ratio_ridge), fill = mle_exists)) +
    geom_text(aes(kappa, sqrt(ratio_ridge) * 1.1,
                  label = sprintf("%.2f", sqrt(ratio_ridge))),
              size = 2,
              angle = 90) +
    facet_grid(. ~ gamma_lab, labeller = label_parsed) +
    geom_hline(aes(yintercept = 1), linetype = 3) +
    labs(x = expression(kappa),
         y = expression(min[alpha]~aMSE[MDYPL]^{1/2} / min[lambda]~aMSE[ridge]^{1/2}),
         fill = "MLE exists") +
    scale_y_continuous(transform = scales::transform_log(),
                       breaks = c(1, 2, 4),
                       limits = c(0.7, 4)) +
    theme_minimal() +
    scale_x_continuous(breaks = (1:9) / 10, limits = c(0, 1.)) +
    theme(legend.position = "right", strip.text = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, size = 7))

mDYPL_ML <- ggplot(all_mses) +
    geom_col(aes(kappa, sqrt(ratio_ML), fill = mle_exists)) +
    geom_text(aes(kappa, sqrt(ratio_ML) * 0.9,
                  label = sprintf("%.2f", sqrt(ratio_ML))),
              size = 2,
              angle = 90) +
    facet_grid(. ~ gamma_lab, labeller = label_parsed) +
    geom_hline(aes(yintercept = 1), linetype = 3) +
    labs(x = NULL,
         y = expression(min[alpha]~aMSE[MDYPL]^{1/2} / aMSE[ML]^{1/2}),
         fill = "MLE exists") +
    scale_y_continuous(transform = scales::transform_log(),
                       breaks = c(1, 2, 4),
                       limits = c(0.7, 4)) +
    theme_minimal() +
    scale_x_continuous(breaks = (1:9) / 10, limits = c(0, 1.)) +
    theme(legend.position = "none",
          axis.text.x = element_blank())

pdf(file.path(figures_path, "rlr_mse_scaled_plots.pdf"), width = 8, height = 5.5)
print(mDYPL_ML / mDYPL_ridge)
dev.off()
