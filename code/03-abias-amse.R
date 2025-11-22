supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "abias-amse.rda")
n_cores <- 10

library("dplyr")
library("ggplot2")
library("patchwork")
library("parallel")
library("brglm2")

source(file.path(supp_path, "code/methods/compute-pt.R"))
source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))

if (file.exists(out_file)) {
    load(out_file)
} else {
    ## Get phase transition curve
    ns <- 200000
    set.seed(123)
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    ga <- c(0.001, 0.01, seq(0, 20, length = 200))
    pt <- compute_pt(gamma_grid = ga, ncores = n_cores, XZU = xzu)

    ## set kappa, gamma grids
    kappas <- seq(0.0075, 0.999, by = 0.005)
    gammas <- seq(0.025, 20, by = 0.1)
    n_kappas <- length(kappas)
    n_gammas <- length(gammas)
    se_pars <- c(0.7, 2, 2)

    result <- mclapply(1:n_gammas, function(ig) {
        gamma = gammas[ig]
        consts <- matrix(NA, n_kappas, 4)
        grads <-  matrix(NA, n_kappas, 3)
        colnames(consts) <- c("mu", "b", "sigma", "max_abs_fun")
        for (ik in 1:n_kappas) {
            kappa = kappas[ik]
            if (ik == 1) {
                res <- solve_se(kappa, gamma, 1/(1 + kappa), start = se_pars)
            } else {
                res <- solve_se(kappa, gamma, 1/(1 + kappa), start = consts[ik - 1, 1:3],
                                init_iter = 0)
            }
            consts[ik, 1:3] <- res
            consts[ik, 4] <- max(abs(attr(res, "funcs")))
            grads[ik, ] <- attr(res, "funcs")
            if (ik %% 10 == 0) {
                cat("MDYPL |",
                    "i_gamma: ", ig, " /", n_gammas, " |",
                    "i_kappa: ", ik, " /", n_kappas, " |",
                    "κ =", round(kappa, digits = 2), ", γ =", round(gamma, digits = 2), " |",
                    round(max(abs(grads[1:ik, ])), digits = 12), "\n")
            }
        }
        data.frame(consts, kappa = kappas, gamma)
    }, mc.cores = n_cores)
    result <- do.call("rbind", result)
    save(pt, result, file = out_file)
}

## Check that all converged
all(result$max_abs_fun < 1e-08)

df <- result |> transform(mse = (sigma^2 + (1 - mu)^2 * gamma^2 / kappa) / 1000)
nbins <- 20

plot_bias <- ggplot(df) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = mu),
                        bins = nbins,
                        show.legend = TRUE) +
    geom_point(data = data.frame(x = 0.2, y = sqrt(0.9)), aes(x, y), pch = 5) +
    geom_line(data = pt, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df$kappa), y = range(df$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression(mu["*"]))

br <- exp(seq(log(min(df$mse)), log(max(df$mse)), length = nbins))
plot_mse <- ggplot(df) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = mse),
                        breaks = br,
                        show.legend = TRUE) +
    geom_point(data = data.frame(x = 0.2, y = sqrt(0.9)), aes(x, y), pch = 5) +
    geom_line(data = pt, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df$kappa), y = range(df$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression(sigma["*"]^2 + kappa^{-1} * gamma^2 * (1 - mu["*"])^2 *~~(x * 10^{-3})))


pdf(file.path(figures_path, "abias-amse.pdf"), width = 12, height = 6)
print(plot_bias + plot_mse + plot_layout(axes = "collect"))
dev.off()
