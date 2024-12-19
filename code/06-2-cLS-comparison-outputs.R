supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

library("dplyr")
library("ggplot2")
library("patchwork")

estimates_files <- dir(results_path, "cLS-MDYPL-comparison-estimates", full.names = TRUE)
consts <- NULL
ests <- NULL
for (f in estimates_files) {
    load(f)
    ests <- rbind(ests, estimates)
    try(consts <- rbind(consts, constants), silent = TRUE)
}

statistics_files <- dir(results_path, "cLS-MDYPL-comparison-statistics", full.names = TRUE)
stats <- NULL
for (f in statistics_files) {
    load(f)
    stats <- rbind(stats, z_statistics)
}

load(file.path(results_path, "cLS-MDYPL-comparison-settings.rda"))

methods_ordered <- c("CLS [O]", "CLS", "rescaled MDYPL [O]", "rescaled MDYPL", "MDYPL")
statistics_ordered <- c("CLS [O]", "CLS", "rescaled MDYPL [O]", "rescaled MDYPL")

stats <- stats |>
    mutate(statistic = factor(statistic, levels = statistics_ordered, ordered = TRUE))
ests <- ests |>
    mutate(method = factor(method, levels = methods_ordered, ordered = TRUE))



root_aMSE <- ests |>
    group_by(method, sample, N, kappa) |> summarize(root_aMSEs = sqrt(mean((estimate - truth)^2))) |>
    group_by(method, N, kappa) |> summarize(root_aMSE = mean(root_aMSEs),
                                            sd = sd(root_aMSEs) / sqrt(n()),
                                            R = n())

aBias <- ests |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())

aBias_zero <- ests |> filter(truth == 0) |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())

aBias_nonzero <- ests |> filter(truth != 0) |>
    group_by(method, sample, N, kappa) |> summarize(aBiases = mean(estimate - truth)) |>
    group_by(method, N, kappa) |> summarize(aBias = mean(aBiases),
                                            sd = sd(aBiases) / sqrt(n()),
                                            R = n())

Bias <- ests |>
    group_by(method, parameter, N, kappa) |> summarize(bias = mean(estimate - truth))


cols0 <- hcl.colors(3, alpha = 0.5)
cols1 <- hcl.colors(3, alpha = 1)
cols <- c(cols0[1], cols1[1], cols0[2], cols1[2])


## root aMSE
p_root_aMSE <- ggplot(root_aMSE |> subset(method != "MDYPL")) +
    geom_point(aes(x = factor(N), y = root_aMSE, col = method, fill = method, shape = method),
               position = position_dodge(0.7)) +
    geom_errorbar(aes(x = factor(N), ymin = root_aMSE - 3 * sd, ymax = root_aMSE + 3 * sd, col = method),
                  position = position_dodge(0.7), width = 0) +
    geom_line(data = settings, aes(x = as.numeric(factor(n)), y = root_aMSE), col = "grey", linetype = 2) +
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


## p-value distribution
## Take the first 25 zero parameters
which_par <- 6:30
which_N <- c(400, 800, 2000)
trans <- hcl.colors(1, alpha = 0)
p_pvalue <- ggplot(stats |> filter(parameter %in% which_par, N %in% which_N)) +
    geom_hline(aes(yintercept = 1), col = "grey") +
    geom_histogram(aes(2 * pnorm(-abs(z)), after_stat(density), col = statistic, fill = statistic),
                   breaks = seq(0, 1, by = 0.05)) +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_color_manual(values = c(cols[2], cols[2], cols[4], cols[4])) +
    scale_fill_manual(values = c(trans, cols[1], trans, cols[3])) +
    facet_grid(kappa + N ~ statistic, labeller = label_bquote(atop(kappa == .(kappa), N == .(N)))) +
    theme_minimal() +
    labs(x = expression(2 * Phi(-"|"~z~"|"))) +
    theme(legend.position = "none")

pdf(file.path(figures_path, "cLS-vs-mDYPL-estimation-inference.pdf"), width = 12, height = 6)

((p_root_aMSE / p_aBias_zero / p_aBias_nonzero) +
 plot_layout(axes = "collect", guides = "collect") & theme(legend.position = "top")) |
    p_pvalue

dev.off()


