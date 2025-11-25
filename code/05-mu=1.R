supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "mu=1.rda")
n_cores <- 10

library("dplyr")
library("ggplot2")
library("patchwork")
library("brglm2")
library("nleqslv")
library("progressr")
handlers("cli")
library("future.apply")
plan(multisession, workers = n_cores)

source(file.path(supp_path, "code/methods/compute-pt.R"))


c2u <- function(pars) {
    c(log(pars[1] / (1 - pars[1])), log(pars[2]), log(pars[3]))
}

u2c <- function(pars) {
    c(1 / (1 + exp(-pars[1])), exp(pars[2:3]))
}

## pars = c(logit(alpha), log(b), log(sigma))
obj_mu_one <- function(pars, kappa, gamma, prox_tol = 1e-10) {
    pars <- u2c(pars)
    alpha <- pars[1]
    b <- pars[2]
    sigma <- pars[3]
    se0(1, b, sigma, kappa, gamma, alpha, prox_tol = prox_tol)
}

## start always in alpha, b, sigma scale
alpha_mu_one <- function(start, kappa, gamma, prox_tol = 1e-10) {
    start <- c2u(start)
    res <- nleqslv(start, obj_mu_one, kappa = kappa, gamma = gamma, prox_tol = prox_tol)
    out <- u2c(res$x)
    attr(out, "linf") <- max(abs(res$fvec))
    out
}

alpha_min_mse <- function(start, kappa, gamma, int = c(0, 1), init_iter = 0,
                          prox_tol = 1e-10) {
    g <- function(alpha) {
        pars <- solve_se(kappa, gamma, alpha, start = start, init_iter = init_iter,
                         control = list(ftol = 1e-12), prox_tol = prox_tol)
        pars[3] / pars[1]
    }
    optimise(g, int, tol = 1e-08)$minimum
}

if (file.exists(out_file)) {
    load(out_file)
} else {
    ## Compute phase transition curve
    ## Get phase transition curve
    ns <- 200000
    set.seed(123)
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    ga <- c(0.001, 0.01, seq(0, 20, length = 200))
    pts <- compute_pt(gamma_grid = ga, ncores = n_cores, XZU = xzu)

    ## Compute alpha such that mu = 1
    kappas <- seq(0.05, 0.975, by = 0.0125)
    n_kappas <- length(kappas)
    gammas <- seq(0.05, 10.0, by = 0.25)
    n_gammas <- length(gammas)
    se_pars <- c(0.8, 2, 2)
    df_unbiased <- NULL
    with_progress({
        every <- 20
        pro <- progressor(n_gammas * n_kappas / every)
        for (ig in 1:n_gammas) {
            gamma <- gammas[ig]
            start <- if (ig == 1) se_pars else consts[1, 1:3]
            consts <- matrix(NA, n_kappas, 4)
            colnames(consts) <- c("alpha", "b", "sigma", "linf")
            for (ik in 1:n_kappas) {
                kappa <- kappas[ik]
                if (ik > 1) start <- consts[ik - 1, 1:3] * c(1.0, 1.4, 1.2)
                res <- alpha_mu_one(start, kappa, gamma, prox = 1e-08)
                consts[ik, 1:3] <- res
                consts[ik, 4] <- attr(res, "linf")
                if (ik %% every == 0)
                    pro(paste("gamma =", gamma, "kappa =", kappa, "|",
                              "linf", round(attr(res, "linf"), 12), "\n"))
            }
            df_unbiased <- rbind(df_unbiased,
                                 data.frame(kappa = kappas, gamma = gamma, consts))
        }
    })

    ## Compute alpha such that mse is minimized
    ## Get solutions at alpha_grid[1]
    alpha_grid <- c(seq(0.01, 0.95, length = 60),
                    seq(0.95, 0.99, length = 40))
    n_a <- length(alpha_grid)
    se_pars <- c(0.8, 0.2, 0.2)
    kg <- expand.grid(kappa = kappas, gamma = gammas)
    n_kg <- nrow(kg)
    start_sol <- matrix(NA, n_kg, 4)
    with_progress({
        every <- 20
        pro <- progressor(n_kg)
        for (i in 1:n_kg) {
            kappa <- kg[i, "kappa"]
            gamma <- kg[i, "gamma"]
            alpha <- alpha_grid[1]
            if ((i - 1) %% length(kappas) == 0) {
                start <- if (i == 1) se_pars else start_sol[i - length(kappas), 1:3]
            } else {
                start <- start_sol[i - 1, 1:3]
            }
            res <- solve_se(kappa, gamma, alpha, start = start,
                            init_iter = 0, control = list(ftol = 1e-12))
            start_sol[i, 1:3] <- res
            start_sol[i, 4] <- max(abs(attr(res, "funcs")))
            if (i %% every == 0)
                pro(paste("kappa =", kappa, "gamma =", gamma, "alpha =", alpha,
                          "linf", round(max(start_sol[, 4], na.rm = TRUE), 12)))
        }
        colnames(start_sol) <- c("mu", "b", "sigma", "linf")
    })

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
                alpha <- alpha_grid[ia]
                start <- if (ia == 1) start else consts[ia - 1, 1:3]
                res <- solve_se(kappa, gamma, alpha, start = start, init_iter = 0,
                            control = list(ftol = 1e-12))
                consts[ia, 1:3] <- res
                consts[ia, 4] <- max(abs(attr(res, "funcs")))
                if (ia %% every == 0) {
                    pro(paste("kappa =", kappa, "gamma =", gamma, "alpha =", alpha, "linf",
                              round(max(consts[, 4], na.rm = TRUE), 12)))
                }
            }
            data.frame(kappa = kappa, gamma = gamma, alpha = alpha_grid, consts)
        }, future.seed = TRUE)
    })
    df_mse <- do.call("rbind", results)
    df_mse$rmse <- with(df_mse, sigma / mu)

    ## Refine
    df_min_mse <- df_mse |>
        group_by(kappa, gamma) |>
        slice_min(order_by = rmse, with_ties = FALSE) |>
        ungroup()
    with_progress({
        pro <- progressor(nrow(df_min_mse))
        alpha_exact <- future_lapply(1:nrow(df_min_mse), function(i) {
            cpars <- df_min_mse[i, ]
            kappa <- cpars$kappa
            gamma <- cpars$gamma
            alpha <- cpars$alpha
            mu <- cpars$mu
            b <- cpars$b
            sigma <- cpars$sigma
            int <- alpha + c(-0.05, 0.05)
            int[int < 0] <- 0.001
            int[int > 1] <- 0.999
            out <- alpha_min_mse(c(mu, b, sigma), kappa, gamma, int = int,
                                 init_iter = 0, prox_tol = 1e-09)
            pro(paste(round(c(out, alpha), 3), "\n"))
            out
        }, future.seed = TRUE)
    })
    df_min_mse$alpha <- unlist(alpha_exact)

    save(df_mse, df_min_mse, df_unbiased, pts, file = out_file)
}

nbins <- 20
br <- seq(min(df_unbiased$alpha), 1, length = nbins)
unbiased <- ggplot(df_unbiased) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = alpha),
                        breaks = br,
                        show.legend = TRUE) +
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df_unbiased$kappa), y = range(df_unbiased$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression({alpha == arg~solve}~paste("{", mu["*"] == 1, "}")))

br <- seq(0, 1, length = nbins)
min_mse <- ggplot(df_min_mse) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = alpha),
                        breaks = br,
                        show.legend = TRUE) +
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df_min_mse$kappa), y = range(df_min_mse$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression({alpha == arg~min}~{sigma["*"]^2 / mu["*"]^2}))


pdf(file.path(figures_path, "alpha-unbiased-min-mse.pdf"), width = 12, height = 6)
print(unbiased + min_mse + plot_layout(axes = "collect"))
dev.off()

