supp_path <- "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")
code_path <- file.path(supp_path, "code")
data_path <- file.path(supp_path, "data")

library("ggplot2")
library("patchwork")

image_files <- dir(results_path, "fou\\+kar", full.names = TRUE)

n_image_files <- length(image_files)
all_results <- as.list(rep(0, n_image_files))
for (i in seq_along(image_files)) {
    im <- image_files[i]
    cat(i, "/", n_image_files, ":", im, "\n")
    load(im)
    all_results[[i]] <- results
}
all_results <- do.call("rbind", all_results)

degrees <- unique(all_results$df)
all_results$gammasq <- all_results$gamma^2

all_results0 <- rbind(all_results |> transform(value = plr, statistic = "PLR"),
                      all_results |> transform(value = rescaled_plr, statistic = "rescaled PLR"))
keep <- c("converged", "statistic", "value", "gammasq", "beta0")

all_results0 <- all_results0[keep] |>
    transform(beta0_fac = paste("beta[0] ==", beta0),
              gammasq_fac = paste("gamma^2 ==", gammasq))

ub <- sort(unique(all_results0$beta0))
ug <- sort(unique(all_results0$gammasq))
all_results0 <- all_results0 |>
    transform(beta0_fac = factor(beta0_fac, paste("beta[0] ==", ub), ordered = TRUE),
              gammasq_fac = factor(gammasq_fac, paste("gamma^2 ==", ug), ordered = TRUE))


qqplots <- ggplot(all_results0 |> subset(converged)) +
    geom_qq(aes(sample = value, col = statistic),
            distribution = stats::qchisq,
            dparams = list(df = degrees),
            alpha = 0.5, size = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1),
                linetype = 2, linewidth = 0.5) +
    facet_grid(gammasq_fac ~ beta0_fac, labeller = label_parsed) +
    scale_color_grey(start = 0.6, end = 0.3) +
    theme_minimal() +
    theme(legend.position = "bottom")

pdf(file.path(figures_path, "case-study-qqplots.pdf"), width = 4, height = 4)
print(qqplots)
dev.off()

