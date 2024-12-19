supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")

load(file.path(results_path, "state_evolution_ML.rda"))
load(file.path(results_path, "state_evolution_mDYPL.rda"))
load(file.path(results_path, "state_evolution_ridge.rda"))


kgambs <- kgambs |> mutate(mse = kappa * sigma^2 / mu^2)
kglmbs <- kglmbs |> mutate(mse = kappa * sigma^2 / mu^2)
kgmbs <- kgmbs |> mutate(mse = kappa * sigma^2 / mu^2)

min_mse_mDYPL <- kgambs |>
    group_by(kappa, gamma) |>
    summarize(alpha_opt = alpha[which.min(mse)],
              mse = min(mse)) |>
    mutate(method = "mDYPL")

min_mse_ridge <- kglmbs |>
    group_by(kappa, gamma) |>
    summarize(lambda_opt = lambda[which.min(mse)],
              mse = min(mse)) |>
    mutate(method = "ridge")

min_mse_ML <- kgmbs |>
    group_by(kappa, gamma) |>
    summarize(mse = min(mse)) |>
    mutate(method = "ML")

all_mses <- rbind(min_mse_ML[c("kappa", "gamma", "mse", "method")],
                  min_mse_mDYPL[c("kappa", "gamma", "mse", "method")],
                  min_mse_ridge[c("kappa", "gamma", "mse", "method")])

all_mses <- all_mses |>
    pivot_wider(names_from = method, values_from = mse) |>
    mutate(ratio_ridge = mDYPL / ridge,
           ratio_ML = mDYPL / ML,
           mle_exists = !is.na(ratio_ML)) |>
    data.frame()

all_mses <- all_mses |> mutate(gamma_lab = factor(paste("gamma ==", gamma),
                                                  levels = paste("gamma ==", unique(gamma)),
                                                  ordered = TRUE))

mDYPL_ridge <- ggplot(all_mses) +
    geom_col(aes(kappa, sqrt(ratio_ridge), fill = mle_exists)) +
    geom_text(aes(kappa, sqrt(ratio_ridge) * 1.1,
                  label = sprintf("%.2f", sqrt(ratio_ridge))),
              size = 2,
              angle = 90) +
    facet_grid(. ~ gamma_lab, labeller = label_parsed) +
    geom_hline(aes(yintercept = 1), linetype = 3) +
    labs(x = expression(kappa),
         y = expression(min[alpha]~aMSE[mDYPL]^{1/2} / min[lambda]~aMSE[ridge]^{1/2}),
         fill = "MLE exists") +
    scale_y_continuous(transform = scales::transform_log(),
                       breaks = c(0.5, 1, 2, 4),
                       limits = c(0.5, 4)) +
    theme_minimal() +
    scale_x_continuous(breaks = (1:9) / 10, limits = c(0, 1.)) +
    theme(legend.position = "bottom", strip.text = element_blank(),
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
         y = expression(min[alpha]~aMSE[mDYPL]^{1/2} / aMSE[ML]^{1/2}),
         fill = "MLE exists") +
    scale_y_continuous(transform = scales::transform_log(),
                       breaks = c(0.5, 1, 2, 4),
                       limits = c(0.5, 4)) +
    theme_minimal() +
    scale_x_continuous(breaks = (1:9) / 10, limits = c(0, 1.)) +
    theme(legend.position = "none",
          axis.text.x = element_blank())


pdf(file.path(figures_path, "rlr_mse_scaled_plots.pdf"), width = 8, height = 6)
print(mDYPL_ML / mDYPL_ridge)
dev.off()


