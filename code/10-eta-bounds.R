supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "upsilon-gamma.rda")
n_cores <- 10

library("parallel")
library("ggplot2")
## library("brglm2")
devtools::load_all("~/Repositories/brglm2")

get_upsilon_gamma <- function(kappa, gammas, alphas, start) {
    n_gamma <- length(gammas)
    upsilon <- rep(NA, n_gamma)
    constants <- matrix(NA, n_gamma, 4)
    colnames(constants) <- c("mu", "b", "sigma", "linf")
    for (j in 1:n_gamma) {
        gamma <- gammas[j]
        alpha <- alphas[j]
        consts <- start <- solve_se(kappa, gamma, alpha, start = start, init_iter = 0)
        constants[j, 1:3] <- consts
        constants[j, 4] <- max(abs(attr(consts, "funcs")))
        upsilon[j] <- sqrt(consts[1]^2 * gamma^2 + kappa * consts[3]^2)
        if (j %% 10 == 0) cat(j, "/", n_gamma, "kappa =", kappa, "gamma =", gamma, "alpha =", alpha, ":", constants[j, ], "\n")
    }
    data.frame(kappa = kappa, gamma = gammas, alpha = alphas, upsilon = upsilon, constants)
}

if (file.exists(out_file)) {
    load(out_file)
} else {
    gammas <- seq(0.001, 15, by = 0.1)
    kappas <- seq(0.1, 0.9, 0.1)
    alphas <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
    n_gammas <- length(gammas)
    n_alphas <- length(alphas)
    n_kappas <- length(kappas)
    ak <- expand.grid(alpha = alphas,
                      kappa = kappas)

    ## Get starting values
    c_results <- data.frame(mu = 0.8, b = 1, sigma = 1)
    starting_values <- NULL
    for (ik in 1:n_kappas) {
        for (ia in 1:n_alphas) {
            kappa <- kappas[ik]
            alpha <- alphas[ia]
            if (ik > 1 & ia == 1) {
                v <- subset(starting_values, kappa == kappas[ik - 1] & alpha == alphas[1])
                start <- unlist(v[1, c("mu", "b", "sigma")]) * 1.2
            } else {
                start <- unlist(c_results[1, c("mu", "b", "sigma")])
            }
            c_results <- get_upsilon_gamma(kappa, gammas[1], alpha, start = start)
            starting_values <- rbind(starting_values, c_results)
        }
    }

    results <- mclapply(1:nrow(ak), function(j) {
        kappa <- ck <- ak[j, "kappa"]
        alpha <- ca <- ak[j, "alpha"]
        start <- unlist(subset(starting_values, kappa == ck & alpha == ca)[1, c("mu", "b", "sigma")])
        get_upsilon_gamma(kappa, gammas, rep(alpha, n_gammas), start = start)
    }, mc.cores = 10)
    results <- do.call("rbind", results)

    ## Adaptive alpha
    alpha_adapt <- plogis(gammas/2)
    start <- c(0.8, 1, 1)
    results_adapt <- mclapply(1:n_kappas, function(ka) {
        kappa <- kappas[ka]
        get_upsilon_gamma(kappa, gammas, alpha_adapt, start = start)
    }, mc.cores = 10)
    results_adapt <- do.call("rbind", results_adapt)
    save(results, results_adapt, file = out_file)
}

pdf(file.path(figures_path, "upsilon-gamma.pdf"), width = 8, height = 2)
ggplot(results) +
   geom_line(aes(gamma, upsilon, color = alpha, group = alpha)) +
   geom_line(data = results_adapt, aes(gamma, upsilon), col = "grey", linewidth = 1) +
   facet_grid(~ kappa, labeller = label_bquote(cols = kappa == .(kappa))) +
   labs(x = expression(gamma), y = expression(upsilon)) +
   scale_color_continuous(type = "viridis") +
   coord_cartesian(y = range(results$upsilon)) +
   theme_minimal()
dev.off()

paths <- dir(results_path, "cLS-MDYPL-comparison-estimates-", full.names = TRUE)
all_estimates <- all_constants <- NULL
for (path in paths) {
    load(path)
    all_estimates <- rbind(all_estimates, estimates)
    all_constants <- rbind(all_constants, constants)
}

bounds <- results |>
    subset(kappa %in% c(0.2, 0.5) & alpha == 0.95) |>
    aggregate(upsilon ~ kappa + alpha, FUN = \(x) c(min = min(x), max = max(x))) |>
    merge(unique(all_constants[c("kappa", "N")]), by = "kappa")
bounds$minimum <- bounds$upsilon[, "min"]
bounds$maximum <- bounds$upsilon[, "max"]



pdf(file.path(figures_path, "cLS-MDYPL-upsilon.pdf"), width = 8, height = 3)
ggplot(all_constants) +
    geom_histogram(aes(upsilon), bins = 50) +
   geom_rect(data = bounds, aes(xmin = minimum, xmax = maximum, ymin = 0, ymax = Inf), alpha = 0.2) +
    facet_grid(kappa ~ N, labeller = label_bquote(cols = N == .(N), rows = kappa == .(kappa))) +
    labs(x = expression(upsilon)) +
    theme_minimal()
dev.off()
