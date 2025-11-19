supp_path <- "."
figures_path <- file.path(supp_path, "figures")
results_path <- file.path(supp_path, "results/new")
out_file <- file.path(results_path, "mfeat-simu.rda")
n_cores <- 10

## library("brglm2")
devtools::load_all("~/Repositories/brglm2")
library("ggplot2")
library("patchwork")
library("parallel")

data("MultipleFeatures", package = "brglm2")


simu_y <- function(X, beta) {
    n <- nrow(X)
    eta <- drop(X %*% beta)
    rbinom(n, 1, plogis(eta))
}

if (file.exists(out_file)) {
    load(out_file)
} else {
    ## Define full and nested models
    full_terms <- c("fou", "kar")
    nest_terms <- c("fou")
    ## Combinations of gamma and beta0 to attempt
    gb0 <- expand.grid(gamma = sqrt(c(1, 2, 4, 8, 16)),
                       beta0 = c(-3, -2, -1, 0))
    vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
    train_id <- which(MultipleFeatures$training)
    MultipleFeatures[train_id, vars] <- scale(MultipleFeatures[train_id, vars], scale = FALSE)
    MultipleFeatures[-train_id, vars] <- scale(MultipleFeatures[-train_id, vars], scale = FALSE)
    MultipleFeatures <- MultipleFeatures |> transform(y = digit == 7)
    full_fm <- formula(paste("y ~", paste(vars, collapse = " + ")))
    nest_vars <- grep("fou", vars, value = TRUE)
    nest_fm <- formula(paste("y ~", paste(nest_vars, collapse = " + ")))
    X_full <- model.matrix(full_fm, data = subset(MultipleFeatures, training))
    X_nest <- model.matrix(nest_fm, data = subset(MultipleFeatures, training))
    full_p <- ncol(X_full) - 1
    nest_p <- ncol(X_nest) - 1
    ## Characteristics
    kappa <- full_p / nrow(X_full)
    alpha <- 1 / (1 + kappa)
    degrees <- full_p - nest_p
    se_start <- expand.grid(mu = c(0.3, 0.5),
                            b = c(2., 1., 0.5),
                            sigma = c(2, 2, 3.),
                            beta0 = c(-3, -2, 0, 2)) |> as.matrix()
    n_start <- nrow(se_start)
    n_simu <- 1000
    results <- NULL
    for (j in seq.int(nrow(gb0))) {
        gamma <- gb0[j, "gamma"]
        beta0 <- gb0[j, "beta0"]
        cat(j, "/", nrow(gb0), ": Computing for gamma^2 =", gamma^2, "beta0 =", beta0, "\n")
        ## Truth
        set.seed(123)
        full_coefs <- structure(rep(0, ncol(X_full)), names = colnames(X_full))
        full_coefs[nest_vars] <- rnorm(nest_p)
        full_coefs[] <- c(beta0, gamma * full_coefs[-1] / sd(X_full[, -1] %*% full_coefs[-1]))
        ## Simulate data sets
        simu_data <- replicate(n_simu, simu_y(X_full, full_coefs), simplify = FALSE)
        results0 <- mclapply(seq.int(n_simu), function(k) {
            df <- data.frame(y = simu_data[[k]], X_full)
            m_full <- glm(full_fm, family = binomial(), data = df, method = mdyplFit, alpha = alpha)
            m_nest <- glm(nest_fm, family = binomial(), data = df, method = mdyplFit, alpha = alpha)
            plr <- plrtest(m_nest, m_full, hd_correction = FALSE)[2, "Deviance"]
            for (att in 1:n_start) {
                rescaled_plr <- try(plrtest(m_nest, m_full, hd_correction = TRUE,
                                            solve_se_control = list(init_iter = 10,
                                                                    start = se_start[att, ],
                                                                    control = list(ftol = 1e-12, xtol = 1e-12))),
                                    silent = TRUE)
                if (inherits(rescaled_plr, "try-error")) {
                    linf <- rescaled_plr <- NA
                    se_pars <- rep(NA, 4)
                    next
                } else {
                    se_pars <- attr(rescaled_plr, "se_parameters")
                    linf <- max(abs(attr(se_pars, "funcs")))
                    if (linf > 1e-07) {
                        rescaled_plr <- NA
                        next
                    } else {
                        rescaled_plr <- rescaled_plr[2, "Deviance"]
                        break
                    }
                }
            }
            msg <- paste0(k, "/", n_simu, " : ", "gamma^2 = ", gamma^2, " beta0 = ", beta0)
            if (att == n_start) cat(msg, ": Failed\n") else cat(msg, ": Done\n")
            data.frame(value = c(plr, rescaled_plr),
                       statistic = c("PLR", "rescaled PLR"),
                       sample = k,
                       mu = se_pars[1],
                       b = se_pars[2],
                       sigma = se_pars[3],
                       theta0hat = se_pars[4],
                       linf = linf,
                       beta0 = beta0,
                       gamma = gamma)
        }, mc.cores = n_cores)
        results <- rbind(results,
                         do.call("rbind", results0))
    }

    beta0_levs <- paste("beta[0] ==", sort(unique(results$beta0)))
    gammasq_levs <- paste("gamma^2 ==", sort(unique(results$gamma^2)))
    results <- results |>
        transform(beta0_fac = factor(paste("beta[0] ==", beta0),
                                     levels = beta0_levs,
                                     ordered = TRUE),
                  gammasq_fac = factor(paste("gamma^2 ==", gamma^2),
                                       levels = gammasq_levs,
                                       ordered = TRUE))

    save(degrees, results, file = out_file)
}


qqplots <- ggplot(subset(results, linf < 1e-07)) +
    geom_qq(aes(sample = value, col = statistic),
            distribution = stats::qchisq,
            dparams = list(df = degrees),
            alpha = 0.5, size = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1),
                linetype = 2, linewidth = 0.5) +
    facet_grid(gammasq_fac ~ beta0_fac, labeller = label_parsed) +
    scale_color_grey(start = 0.6, end = 0.3) +
    theme_minimal() +
    theme(legend.position = "top")

pdf(file.path(figures_path, "case-study-qqplots.pdf"), width = 4, height = 4)
print(qqplots)
dev.off()
