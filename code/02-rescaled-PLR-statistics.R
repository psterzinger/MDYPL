supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results")
out_file <- file.path(results_path, "qqplots-rescaled-plr.rda")
n_cores <- 10

library("dplyr")
library("ggplot2")
library("patchwork")
library("parallel")
library("RcppNumerical")
## library("brglm2")
devtools::load_all("~/Repositories/brglm2")

source(file.path(supp_path, "code/methods/compute-pt.R"))
source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))
source(file.path(supp_path, "code/methods/fit-mdypl.R"))

if (file.exists(out_file)) {
    load(out_file)
} else {
    gamma <- sqrt(5)
    ## Get phase transition curve
    ns <- 200000
    set.seed(123)
    xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
    ga <- c(0.001, 0.01, seq(0, 20, length = 200))
    pt_point <- compute_pt(gamma_grid = gamma, ncores = n_cores, XZU = xzu)

    kga <- expand.grid(kappa = c(0.1, 0.5),
                       gamma = gamma,
                       alpha = c(1, 3/4, 1/2, 1/4))
    kga <- rbind(kga,
                 c(0.1, sqrt(5), 1 / (1 + 0.1)),
                 c(0.5, sqrt(5), 1 / (1 + 0.5)))
    mle_exists <- kga$kappa < pt_point$kappa
    mbs <- matrix(NA, ncol = 3, nrow = nrow(kga))
    se_pars <- c(0.7, 2, 2)
    for (s in 1:nrow(kga)) {
        kappa <- kga[s, "kappa"]
        gamma <- kga[s, "gamma"]
        alpha <- kga[s, "alpha"]
        if (alpha == 1 & !mle_exists[s]) next
        mbs[s, ] <- solve_se(kappa, gamma, alpha, start = se_pars, init_iter = 10)
        cat("kappa =", kappa, "gamma =", gamma, "alpha", alpha, "Done.\n")
    }
    colnames(mbs) <- c("mu", "b", "sigma")
    kgambs <- cbind(kga, mbs)


    R <- 1000
    n <- 2000

    ## Compute PLR and rescaled PLR statistics
    set.seed(123)
    results <- NULL
    for (s in 1:nrow(kgambs)) {
        kappa <- kgambs[s, "kappa"]
        alpha <- kgambs[s, "alpha"]
        cat("kappa =", kappa, "gamma =", gamma, "alpha", alpha, "\n")
        ## move on if ML (alpha = 1) and MLE does not exist
        if (alpha == 1 & !mle_exists[s]) next
        mu <- kgambs[s, "mu"]
        b <- kgambs[s, "b"]
        sigma <- kgambs[s, "sigma"]
        p <- n * kappa
        beta0 <- rep(c(0, 1), each = p / 2)
        ## resclale so that signal strength is gamma
        beta0 <- sqrt(n) * gamma * beta0 / sqrt(sum(beta0^2))
        seeds <- generate_unique_seeds(R)
        res <- mclapply(1:R, function(i) {
            if (i %% 100 == 0) cat(i, "out of", R, "\n")
            set.seed(seeds[i])
            X <- matrix(rnorm(n * p), nrow = n, ncol = p) / sqrt(n)
            probs <- plogis(c(X %*% beta0))
            y <- rbinom(n, 1, probs)
            dat <- data.frame(y, X)
            ## fm_full <- paste0("y ~ -1 +", paste0("X", 1:p, collapse = " + "))
            ## fm_nest1 <- paste0("y ~ -1 +", paste0("X", 6:p, collapse = " + "))
            ## fm_nest2 <- paste0("y ~ -1 +", paste0("X", 51:p, collapse = " + "))
            ## m_full <- glm(fm_full, family = binomial(), data = dat, method = "mdypl_fit", alpha = alpha,
            ##               start = beta0)
            ## m_nest1 <- glm(fm_nest1, family = binomial(), data = dat, method = "mdypl_fit", alpha = alpha,
            ##                start = beta0[-c(1:5)])
            ## m_nest2 <- glm(fm_nest2, family = binomial(), data = dat, method = "mdypl_fit", alpha = alpha,
            ##                start = beta0[-c(1:50)])
            m_full <- fit_mdypl(X, y, alpha = alpha, start = beta0)
            m_nest1 <- fit_mdypl(X[, -c(1:5)], y, alpha = alpha, start = beta0[-c(1:5)])
            m_nest2 <- fit_mdypl(X[, -c(1:50)], y, alpha = alpha, start = beta0[-c(1:50)])
            stat1 <- 2 * (m_full$pl - m_nest1$pl)
            stat2 <- 2 * (m_full$pl - m_nest2$pl)
            r_stat1 <- stat1 * b / (kappa * sigma^2)
            r_stat2 <- stat2 * b / (kappa * sigma^2)
            data.frame(value = c(stat1, stat2, r_stat1, r_stat2),
                       method = rep(c("PLR", "rescaled_PLR"), each = 2),
                       kappa = kappa,
                       gamma = gamma,
                       alpha = alpha,
                       df = rep(c(5, 50), 2))
        }, mc.cores = n_cores)
        results <- rbind(results, do.call("rbind", res))
    }
    results <- results |>
        transform(alpha_fac = as.character(MASS::fractions(alpha)),
                  kappa_fac = paste("kappa ==", kappa))
    results[results$alpha_fac %in% c("10/11", "2/3"), "alpha_fac"] <- "1 / (1 + kappa)"
    results$alpha_fac <- factor(results$alpha_fac, levels = c("1", "3/4", "1/2", "1/4", "1 / (1 + kappa)"), ordered = TRUE)
    levels(results$alpha_fac) <- paste("alpha ==", levels(results$alpha_fac))
    save(kgambs, results, file = out_file)
}

qq_5 <- ggplot(results |> subset(df == 5)) +
    geom_qq(aes(sample = value, col = method), distribution = qchisq, dparams = list(df = 5),
           alpha = 0.5, size = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    facet_grid(alpha_fac ~ kappa_fac, labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_grey(start = 0.6, end = 0.3, name = "statistic") +
    labs(title = expression({{beta[paste(0, ",", 1)] == ldots} == beta[paste(0, ",", 5)]} == 0))

qq_50 <- ggplot(results |> subset(df == 50)) +
   geom_qq(aes(sample = value, col = method), distribution = qchisq, dparams = list(df = 50),
           alpha = 0.5, size = 0.5) +
   geom_abline(aes(intercept = 0, slope = 1)) +
   facet_grid(alpha_fac ~ kappa_fac, labeller = label_parsed) +
   theme_minimal() +
   theme(legend.position = "right") +
   scale_color_grey(start = 0.6, end = 0.3, name = "statistic") +
   labs(title = expression({{beta[paste(0, ",", 1)] == ldots} == beta[paste(0, ",", 50)]} == 0))

pdf(file.path(figures_path, "qqplots-rescaled-plr.pdf"), width = 10, height = 5)
print(qq_5 + qq_50 + plot_layout(axes = "collect"))
dev.off()
