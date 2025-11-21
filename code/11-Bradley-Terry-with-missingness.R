supp_path = "."
figures_path <- file.path(supp_path, "figures")
results_path = file.path(supp_path, "results/new")
n_cores <- 10

library("parallel")
library("ggplot2")
## library("brglm2")
devtools::load_all("~/Repositories/brglm2")
library("detectseparation")

source(file.path(supp_path, "code/methods/generate-unique-seeds.R"))

simulate_BT_df <- function(abilities, gamma = 1, kappa = 0.1) {
    n_players <- length(abilities)
    players <- names(abilities)
    if (is.null(players)) {
        width <- nchar(as.character(n_players))
        players <- sprintf(paste0("P%0", width, "d"), 1:n_players)
    }
    player1 <- rep(players, each = n_players)
    player2 <- rep(players, n_players)
    ids <- player1 != player2
    competitions <- data.frame(player1, player2) |> subset(ids)
    X <- matrix(0, nrow(competitions), n_players)
    colnames(X) <- players
    for (j in seq.int(nrow(X))) {
        player1 <- competitions[j, "player1"]
        player2 <- competitions[j, "player2"]
        X[j, player1] <- -1
        X[j, player2] <- 1
    }
    abilities <- gamma * abilities / drop(sd(X %*% abilities))
    probs <- plogis(outer(abilities, abilities, "-"))
    rownames(probs) <- colnames(probs) <- players
    outcome <- rbinom(n_players^2, 1, c(probs))[ids]
    out <- data.frame(y = outcome, X)
    keep_p <- 1 / (n_players * kappa)
    keep <- sample(1:nrow(out), ceiling(nrow(out) * keep_p),
                   replace = FALSE)
    out <- out[keep, ]
    attr(out, "abilities") <- abilities
    attr(out, "keep") <- keep_p
    class(out) <- c("BT_df", class(out))
    out
}

n_simu <- 200
kappa0 <- 0.2
gamma0 <- 4
alpha0 <- 0.7

out_file <- file.path(results_path, paste0("BT-MAR-kappa", kappa0, "-gamma-",
                                           gamma0, "-alpha-", alpha0, "-R-",
                                           n_simu, ".rda"))

if (file.exists(out_file)) {
    load(out_file)
} else {
    abilities0 <- c(rep(0, 10), seq(-100, 100, length = 290)) / 20
    n_players <- length(abilities0)
    ## Set up formulae
    width <- nchar(as.character(n_players))
    names(abilities0) <- sprintf(paste0("P%0", width, "d"), 1:n_players)
    fixed <- names(which(abilities0 == 0))
    i_fixed <- names(abilities0) == fixed[1]
    BT_formula <- formula(paste("y ~ -1 + . -", fixed[1]))
    BT_formula_nest <- formula(paste("y ~ -1 + . -", paste(fixed, collapse = "-")))
    set.seed(123)
    seeds <- generate_unique_seeds(n_simu)
    res <- mclapply(1:n_simu, function(j) {
        set.seed(seeds[j])
        dat <- simulate_BT_df(abilities0, kappa = kappa0, gamma = gamma0)
        fit_mdypl <- glm(BT_formula, data = dat, family = binomial(),
                         method = "mdypl_fit",
                         alpha = alpha0)
        is_sep <- glm(BT_formula, data = dat, family = binomial(),
                      method = "detect_separation")
        summ_hd <- try(summary(fit_mdypl, hd_correction = TRUE), silent = TRUE)
        se_failed <- inherits(summ_hd, "try-error")
        coefs_mdypl <-  na.omit(coef(fit_mdypl))
        fit_mdypl_nest <- glm(BT_formula_nest, data = dat, family = binomial(),
                              method = "mdypl_fit",
                              alpha = alpha0)
        pv <- plrtest(fit_mdypl_nest, fit_mdypl)[2, "Pr(>Chi)"]
        if (se_failed) {
            ss <- NA
            coefs_sc_mdypl <- rep(NA, length(coefs_mdypl))
            se_pars <- rep(NA, 3)
            sc_pv <- NA
        } else {
            ss <- summ_hd$signal_strength
            coefs_sc_mdypl <- coef(summ_hd)[, "Estimate"]
            se_pars <- summ_hd$se_parameters
            sc_pv <- plrtest(fit_mdypl_nest, fit_mdypl, hd_correction = TRUE,
                             se_start = se_pars)[2, "Pr(>Chi)"]
        }
        cat(j, ":", if (se_failed) "failed" else "done", "ss", ":", ss, "\n")
        co <- co_sc <- is_inf <- rep(NA, n_players)
        names(co) <- names(co_sc) <- names(is_inf) <- colnames(dat)[-1]
        co[names(coefs_mdypl)] <- coefs_mdypl
        co_sc[names(coefs_mdypl)] <- coefs_sc_mdypl
        is_inf[names(coef(is_sep))] <- is.infinite(coef(is_sep))
        list(coefs_mdypl = co,
             coefs_sc_mdypl = co_sc,
             inf_coefs = is_inf,
             se_pars = se_pars,
             pv = pv,
             sc_pv = sc_pv,
             ss = ss,
             separation = is_sep$outcome)
    },
    mc.cores = n_cores)
    abilities <- attr(simulate_BT_df(abilities0, kappa = kappa0, gamma = gamma0), "abilities")
    save(n_players, res, abilities, file = out_file)
}

is_separated <- sapply(res, "[[", "separation")
inf_coefs <- sapply(res, "[[", "inf_coefs")
hd_converged <- sapply(res, function(x) !(is.na(x$ss) || max(abs(attr(x$se_pars, "funcs"))) > 1e-08))

## All converged?
all(hd_converged)

## Number separated
sum(is_separated)

inf_coefs <- sapply(res[hd_converged], "[[", "inf_coefs")
coefs_mdypl <- sapply(res[hd_converged], "[[", "coefs_mdypl")
coefs_sc_mdypl <- sapply(res[hd_converged], "[[", "coefs_sc_mdypl")
pv <- sapply(res[hd_converged], "[[", "pv")
sc_pv <- sapply(res[hd_converged], "[[", "sc_pv")

## Bias
e_mdypl <- rowMeans(coefs_mdypl, na.rm = TRUE)
e_sc_mdypl <- rowMeans(coefs_sc_mdypl, na.rm = TRUE)
bias_mdypl <- e_mdypl - abilities[names(e_mdypl)]
bias_sc_mdypl <- e_sc_mdypl - abilities[names(e_sc_mdypl)]

## Centered estimates
c_estimates <- data.frame(estimate = c(coefs_mdypl, coefs_sc_mdypl),
                          method = rep(c("MDYPL", "rescaled MDYPL"), each = prod(dim(coefs_mdypl))),
                          ability = c(rep(abilities, ncol(coefs_mdypl)), rep(abilities, ncol(coefs_sc_mdypl))),
                          player = rep(1:nrow(coefs_mdypl), ncol(coefs_mdypl) + ncol(coefs_sc_mdypl)))
biases <- data.frame(bias = c(bias_mdypl, bias_sc_mdypl),
                     method = rep(c("MDYPL", "rescaled MDYPL"), each = n_players),
                     ability = rep(abilities, 2),
                     player = rep(1:nrow(coefs_mdypl), 2))

ests <- ggplot(c_estimates) +
    geom_boxplot(aes(player, estimate - ability, group = player, color = ability),
                 outlier.alpha = 0.2) +
    geom_point(data = biases, aes(player, bias), alpha = 0.1) +
    geom_abline(aes(intercept = 0, slope = 0), col = "grey") +
    facet_grid(method ~ .) +
    labs(x = "player (j)", y = expression(hat(beta)[j]^{DY} - beta[paste(0, ",", j)])) +
    theme_minimal() +
    scale_color_continuous(type = "viridis",
                           name = expression(paste("ability (", beta[paste(0, ",", j)], ")")))

pdf(file.path(figures_path, "bt-centred-estimates.pdf"), width = 9, height = 4)
print(ests)
dev.off()

## Infinite ML estimates
inf_ml_estimates <- data.frame(player = 1:nrow(inf_coefs),
                               estimate = c(inf_coefs),
                               ability = rep(abilities, ncol(inf_coefs))) |>
    aggregate(estimate ~ player + ability, mean)

inf_ests <- ggplot(inf_ml_estimates) +
    geom_col(aes(x = player, y = estimate, col = ability)) +
    theme_minimal() +
    labs(x = "player (j)", y = expression(paste("% of infinite values for ", hat(beta)[j]))) +
    lims(y = c(0, 1)) +
    scale_color_continuous(type = "viridis", name = expression(paste("ability (", beta[paste(0, ",", j)], ")")))

pdf(file.path(figures_path, "bt-inf.pdf"), width = 9, height = 2)
print(inf_ests)
dev.off()


## Equality of abilities
degrees <- sum(abilities == 0) - 1
results <- data.frame(value = c(qchisq(pv, degrees), qchisq(sc_pv, degrees)),
                      statistic = rep(factor(c("PLR", "adjusted PLR"), levels = c("PLR", "adjusted PLR"), ordered = TRUE), each = length(pv)))

qqplots <- ggplot(results) +
    geom_qq(aes(sample = value, col = statistic),
            distribution = stats::qchisq,
            dparams = list(df = degrees),
            alpha = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1),
                linetype = 2, linewidth = 0.5) +
    scale_color_grey(start = 0.6, end = 0.3) +
    lims(x = c(0, 40), y = c(0, 40)) +
   labs(title = expression({{beta[paste(0, ",", 2)] == ldots} == beta[paste(0, ",", 10)]} == 0))  +
    facet_grid(~ statistic) +
    theme_minimal()

pdf(file.path(figures_path, "bt-qqplots.pdf"), width = 9, height = 3)
qqplots
dev.off()
