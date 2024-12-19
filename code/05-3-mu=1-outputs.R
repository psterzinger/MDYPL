supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

library("ggplot2")
library("patchwork")

nbins <- 20
load(file.path(results_path, "mu=1.rda"))
br <- seq(min(df$alpha), 1, length = nbins)
alpha_unbiased <- ggplot(df) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = alpha),
                        breaks = br,
                        show.legend = TRUE) +
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df$kappa), y = range(df$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression({alpha == arg~solve}~paste("{", mu["*"] == 1, "}")))

load(file.path(results_path, "alpha-min-mse.rda"))
br <- seq(0, 1, length = nbins)
alpha_min_mse <- ggplot(df) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = alpha),
                        breaks = br,
                        show.legend = TRUE) +
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df$kappa), y = range(df$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression({alpha == arg~min}~{sigma["*"]^2 / mu["*"]^2}))


pdf(file.path(figures_path, "alpha-unbiased-min-mse.pdf"), width = 12, height = 6)
print(alpha_unbiased + alpha_min_mse + plot_layout(axes = "collect"))
dev.off()

