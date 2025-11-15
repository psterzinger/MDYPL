supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results/new")
n_cores <- 10

## library("brglm2")
devtools::load_all("~/Repositories/brglm2")
library("mvtnorm")
library("dplyr")
library("ggplot2")
library("patchwork")
library("parallel")
library("glmnet")
library("tictoc")


source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))

scaled_beta <- function(n, kappa, gamma) {
    p <- floor(n * kappa)
    Sigma <- (matrix(1, p, p) + diag(p)) / 2
    beta <- rep(0, p)
    beta[1:5] <- 1
    gamma * beta / sqrt(drop(t(beta) %*% Sigma %*% beta))
}

eta_lasso <- function(X, y) {
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)
    best_lambda <- cv_fit$lambda.min
    best_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)
    coefs <- coef(best_model)
    predict(best_model, newx = X, type = "link")
}

simu_X <- function(n, kappa) {
    p <- floor(n * kappa)
    Sigma <- (matrix(1, p, p) + diag(p)) / 2
    rmvnorm(n, rep(0, p), sigma = Sigma)
}

simu_y <- function(X, beta) {
    n <- nrow(X)
    eta <- drop(X %*% beta)
    rbinom(n, 1, exp(eta) / (1 + exp(eta)))
}

cls <- function(X, y, eta, qrX = NULL) {
    eta_norm <- sum(eta^2)
    eta_tanh <- tanh(eta / 2)
    c_hat <- sum(eta * eta_tanh) / eta_norm
    if (is.null(qrX)) {
        qrX <- qr(X)
    }
    beta_hat <- solve.qr(qrX, 2 * y - 1)
    r1 <- X %*% solve.qr(qrX, eta_tanh)
    r2 <- eta %*% crossprod(eta, eta_tanh) / eta_norm
    delta_hat <- solve.qr(qrX, r1 - r2)
    (beta_hat - delta_hat) / c_hat
}

cls_z_stat <- function(estimate, beta0, X, eta, qrX = NULL) {
    eta_tanh <- tanh(eta / 2)
    sigma_hat <- sum(eta * eta_tanh) / sum(eta^2)
    w <- sqrt(1 - eta_tanh^2)    # length n
    if (is.null(qrX)) {
        qrX <- qr(X)
    }
    dX <- backsolve(qr.R(qrX), t(w * qr.Q(qrX)))
    denom <- sqrt(rowSums(dX^2))
    sigma_hat * (estimate - beta0) / denom
}

mdypl <- function(X, y, alpha, start = NULL, se_start = c(0.5, 1, 1)) {
    fit <- glm(y ~ -1 + X, family = binomial(), method = "mdypl_fit", alpha = alpha, start = start)
    ss <- try(summary(fit, hd_correction = TRUE,  se_start = se_start, init_iter = 0), silent = TRUE)
    if (inherits(ss, "try-error")) { # try another time
        ss <- summary(fit, hd_correction = TRUE,  se_start = se_start * 0.9, init_iter = 0)
    }
    consts <- ss$se_parameters
    list(estimate = coef(ss)[, "Estimate"],
         unscaled_estimate = coef(ss)[, "Estimate"] * consts[1],
         z_statistic = coef(ss)[, "z value"],
         consts = consts,
         upsilon = sqrt(ss$signal_strength))
}

R <- 1000
gamma <- 3
alpha <- 0.95
ns <- c(400, 800, 1200, 1600, 2000)
kappas <- c(0.2, 0.5)
nk <- expand.grid(n = ns, kappa = kappas)
nk <- nk |>
    mutate(path = file.path(results_path,
                            paste0("cLS-MDYPL-comparison-estimates-", n, "-", kappa, ".rda")))
nk$mu <- nk$b <- nk$sigma <- NA
n_settings <- nrow(nk)

if (all(file.exists(nk$path))) {
    ## Gather results form images
    all_estimates <- all_z_statistics <- all_constants <- NULL
    for (path in nk$path) {
        load(path)
        all_estimates <- rbind(all_estimates, estimates)
        all_z_statistics <- rbind(all_z_statistics, z_statistics)
        all_constants <- rbind(all_constants, constants)
    }
} else {
    ## General unique seeds
    set.seed(123)
    seeds <- generate_unique_seeds(R * n_settings)
    seeds <- matrix(seeds, ncol = R)
    methds <- c("CLS", "CLS [O]", "rescaled MDYPL", "rescaled MDYPL [O]")
    n_methds <- length(methds)
    stats <-  c("CLS", "CLS [O]", "rescaled MDYPL", "rescaled MDYPL [O]")
    n_stats <- length(stats)
    for (i in 1:n_settings) {
        n <- nk[i, "n"]
        kappa <- nk[i, "kappa"]
        ## Get oracle constants and root MSE
        se_init <- c(0.9, 2, 2)
        consts <- solve_se(nk$kappa[i], gamma, alpha, start = se_init)
        nk$mu[i] <- mu <- consts[1]
        nk$b[i] <- b <- consts[2]
        nk$sigma[i] <- sigma <- consts[3]
        nk$root_aMSE[i] <- sqrt(2 / (n + 1 / kappa)) * sigma / mu
        cat("\nSetting", i, ": n =", n, " kappa=", kappa, "\n")
        set.seed(100)
        X <- simu_X(n, kappa)
        beta_oracle <- scaled_beta(n, kappa, gamma)
        p <- length(beta_oracle)
        eta_oracle <- drop(X %*% beta_oracle)
        ## responses
        y <- mclapply(1:R, function(r) { set.seed(seeds[i, r]); simu_y(X, beta_oracle) },
                      mc.cores = n_cores)
        y <- do.call("cbind", y)
        ## LASSO etas
        tic("Estimating eta")
        eta_hat <- mclapply(1:R, function(r) { set.seed(seeds[i, r]); eta_lasso(X, y[, r]) },
                            mc.cores = n_cores)
        eta_hat <- do.call("cbind", eta_hat)
        toc()
        ## MDYPL
        tic("Computing MDYPL estimates and statistics")
        ## Oracle constants
        consts_oracle <- c(mu, b, sigma)
        tau_oracle <- sqrt((p + 1) / (2 * p))
        ## Get MDYPL fits
        mdypl_fits <- mclapply(1:R, function(r) mdypl(X, y[, r], alpha, se_start = consts_oracle * 0.9),
                               mc.cores = n_cores)
        ## Corrected estimates and z statistics for beta = 0
        beta_mdypl <- do.call("cbind", lapply(mdypl_fits, "[[", "estimate"))
        z_mdypl <- do.call("cbind", lapply(mdypl_fits, "[[", "z_statistic"))
        ## Oracle corrected estimates and z statistics for beta = 0
        beta_mdypl_unscaled <- do.call("cbind", lapply(mdypl_fits, "[[", "unscaled_estimate"))
        beta_mdypl_oracle <- beta_mdypl_unscaled / consts_oracle[1]
        z_mdypl_oracle <- sqrt(n) * tau_oracle * beta_mdypl_unscaled / consts_oracle[3]
        toc()
        ## CLS
        ## Corrected estimates and z statistics for beta = 0
        tic("Compputing CLS estimates and statistics")
        qr_x <- qr(X)
        tic()
        beta_cls <- mclapply(1:R, function(r) cls(X, y[, r], eta_hat[, r], qrX = qr_x),
                             mc.cores = n_cores)
        beta_cls <- do.call("cbind", beta_cls)
        z_cls <- mclapply(1:R, function(r) cls_z_stat(beta_cls[, r], rep(0, p), X, eta_hat[, r], qrX = qr_x),
                          mc.cores = n_cores)
        z_cls <- do.call("cbind", z_cls)
        ## Oracle corrected estimates and z statistics for beta = 0
        beta_cls_oracle <- mclapply(1:R, function(r) cls(X, y[, r], eta_oracle, qrX = qr_x),
                                    mc.cores = n_cores)
        beta_cls_oracle <- do.call("cbind", beta_cls_oracle)
        z_cls_oracle <- mclapply(1:R, function(r) cls_z_stat(beta_cls_oracle[, r], rep(0, p), X, eta_oracle, qrX = qr_x),
                                 mc.cores = n_cores)
        z_cls_oracle <- do.call("cbind", z_cls_oracle)
        toc()
        ## Collect results
        tic("Collecting results")
        ## Constants
        consts <- lapply(mdypl_fits, function(x) {
            co <- x$consts
            data.frame(mu = co[1], b = co[2], sigma = co[3], linf = max(abs(attr(co, "funcs"))),
                       uspilon = x$upsilon, gamma = gamma, kappa = kappa, N = n, p = p)
        })
        constants <- do.call("rbind", consts) |> mutate(sample = 1:R)
        estimates <- data.frame(estimate = c(beta_cls, beta_cls_oracle, beta_mdypl, beta_mdypl_oracle),
                                method = rep(methds, each = R * p),
                                truth = beta_oracle,
                                parameter = rep(1:p, R * n_methds),
                                sample = rep(rep(1:R, each = p), n_methds),
                                kappa = kappa, gamma = gamma, N = n, p = p)
        z_statistics <- data.frame(z = c(z_cls, z_cls_oracle, z_mdypl, z_mdypl_oracle),
                                   statistic = rep(stats, each = R * p),
                                   truth = beta_oracle,
                                   parameter = rep(1:p, R * n_stats),
                                   sample = rep(rep(1:R, each = p), n_stats),
                                   kappa = kappa, gamma = gamma, N = n, p = p)
        toc()
        save(nk, constants, estimates, z_statistics, file = nk$path[i])
    }
}

methods_ordered <- c("CLS [O]", "CLS", "rescaled MDYPL [O]", "rescaled MDYPL", "MDYPL")
statistics_ordered <- c("CLS [O]", "CLS", "rescaled MDYPL [O]", "rescaled MDYPL")
all_z_statistics <- all_z_statistics |>
    mutate(statistic = factor(statistic, levels = statistics_ordered, ordered = TRUE))
all_estimates <- all_estimates |>
    mutate(method = factor(method, levels = methods_ordered, ordered = TRUE))

root_aMSE <- all_estimates |>
    group_by(method, sample, N, kappa) |> summarize(root_aMSEs = sqrt(mean((estimate - truth)^2))) |>
    group_by(method, N, kappa) |> summarize(root_aMSE = mean(root_aMSEs),
                                            sd = sd(root_aMSEs) / sqrt(n()),
                                            R = n())
aBias <- all_estimates |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())
aBias_zero <- all_estimates |> filter(truth == 0) |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())
aBias_nonzero <- all_estimates |> filter(truth != 0) |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())
Bias <- all_estimates |>
    group_by(method, parameter, N, kappa) |> summarize(bias = mean(estimate - truth))

cols0 <- hcl.colors(3, alpha = 0.5)
cols1 <- hcl.colors(3, alpha = 1)
cols <- c(cols0[1], cols1[1], cols0[2], cols1[2])

## root aMSE
p_root_aMSE <- ggplot(root_aMSE) +
    geom_point(aes(x = factor(N), y = root_aMSE, col = method, fill = method, shape = method),
               position = position_dodge(0.7)) +
    geom_errorbar(aes(x = factor(N), ymin = root_aMSE - 3 * sd, ymax = root_aMSE + 3 * sd, col = method),
                  position = position_dodge(0.7), width = 0) +
    geom_line(data = nk, aes(x = as.numeric(factor(n)), y = root_aMSE), col = "grey", linetype = 2) +
    facet_grid(~ kappa, labeller = label_bquote(cols = kappa == .(kappa))) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_shape_manual(values = 20 + c(1, 2, 1, 2)) +
    labs(x = "n", y = "root aMSE") +
    theme_minimal()

## Bias non-zero
p_aBias_nonzero <- ggplot(aBias_nonzero |> subset(method != "MDYPL")) +
    geom_hline(aes(yintercept = 0), col = "grey", linetype = 2) +
    geom_point(aes(x = factor(N), y = aBias, col = method, fill = method, shape = method),
               position = position_dodge(0.7)) +
    geom_errorbar(aes(x = factor(N), ymin = aBias - 3 * sd, ymax = aBias + 3 * sd, col = method),
                  position = position_dodge(0.7), width = 0) +
    facet_grid(~ kappa, labeller = label_bquote(cols = kappa == .(kappa))) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_shape_manual(values = 20 + c(1, 2, 1, 2)) +
    lims(y = c(-0.2, 0.2)) +
    labs(x = "n", y = expression(paste("aBias (", beta[j] != 0, ")"))) +
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())

p_aBias_zero <- ggplot(aBias_zero |> subset(method != "MDYPL")) +
    geom_hline(aes(yintercept = 0), col = "grey", linetype = 2) +
    geom_point(aes(x = factor(N), y = aBias, col = method, fill = method, shape = method),
               position = position_dodge(0.7)) +
    geom_errorbar(aes(x = factor(N), ymin = aBias - 3 * sd, ymax = aBias + 3 * sd, col = method),
                  position = position_dodge(0.7), width = 0) +
    facet_grid(~ kappa, labeller = label_bquote(cols = kappa == .(kappa))) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_shape_manual(values = 20 + c(1, 2, 1, 2)) +
    lims(y = c(-0.02, 0.02)) +
    labs(x = "n", y = expression(paste("aBias (", beta[j] == 0, ")"))) +
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())

trans <- hcl.colors(1, alpha = 0)
which_par <- 6:30
which_N <- c(400, 800, 2000)
p_pvalue <- ggplot(all_z_statistics |> subset(parameter %in% which_par & N %in% which_N)) +
    geom_hline(aes(yintercept = 1), col = "grey") +
    geom_histogram(aes(2 * pnorm(-abs(z)), after_stat(density), col = statistic, fill = statistic),
                   breaks = seq(0, 1, by = 0.05)) +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_color_manual(values = c(cols[2], cols[2], cols[4], cols[4])) +
    scale_fill_manual(values = c(trans, cols[1], trans, cols[3])) +
    lims(y = c(0, 1.3)) +
    facet_grid(kappa + N ~ statistic, labeller = label_bquote(atop(kappa == .(kappa), N == .(N)))) +
    theme_minimal() +
    labs(x = expression(2 * Phi(-"|"~z~"|"))) +
    theme(legend.position = "none")

pdf(file.path(figures_path, "cLS-vs-mDYPL-estimation-inference.pdf"), width = 12, height = 6)
((p_root_aMSE / p_aBias_zero / p_aBias_nonzero) +
 plot_layout(axes = "collect", guides = "collect") & theme(legend.position = "top")) |
    p_pvalue
dev.off()


