supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

load(file.path(results_path, "qqplots-rescaled-plr.rda"))

library("ggplot2")
library("patchwork")

qq_5 <- ggplot(lr_5 |> na.omit()) +
    geom_qq(aes(sample = value, col = method), distribution = qchisq, dparams = list(df = 5),
           alpha = 0.5, size = 0.5) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    facet_grid(alpha_fac ~ kappa_fac, labeller = label_parsed) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_grey(start = 0.6, end = 0.3) +
    labs(title = expression({{beta[paste(0, ",", 1)] == ldots} == beta[paste(0, ",", 5)]} == 0))
qq_50 <- ggplot(lr_50 |> na.omit()) +
   geom_qq(aes(sample = value, col = method), distribution = qchisq, dparams = list(df = 50),
           alpha = 0.5, size = 0.5) +
   geom_abline(aes(intercept = 0, slope = 1)) +
   facet_grid(alpha_fac ~ kappa_fac, labeller = label_parsed) +
   theme_minimal() +
   theme(legend.position = "right") +
   scale_color_grey(start = 0.6, end = 0.3) +
   labs(title = expression({{beta[paste(0, ",", 1)] == ldots} == beta[paste(0, ",", 50)]} == 0))


pdf(file.path(figures_path, "qqplots-rescaled-plr.pdf"), width = 10, height = 5)
print(qq_5 + qq_50 + plot_layout(axes = "collect"))
dev.off()
