supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

load(file.path(results_path, "abias-amse.rda"))

library("ggplot2")
library("patchwork")

nbins <- 20

df <- df |> transform(mse = (sigma^2 + (1 - mu)^2 * gamma^2 / kappa)/1000)

plot_bias <- ggplot(df) +
    geom_contour_filled(aes(x = kappa, y = gamma, z = mu),
                        bins = nbins,
                        show.legend = TRUE) +
    geom_point(data = data.frame(x = 0.2, y = sqrt(0.9)), aes(x, y), pch = 5) +
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
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
    geom_line(data = pts, aes(kappa, gamma), col = "darkgrey", lwd = 1.2) +
    lims(x = range(df$kappa), y = range(df$gamma)) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(x = expression(kappa), y = expression(gamma),
         title = expression(sigma["*"]^2 + kappa^{-1} * gamma^2 * (1 - mu["*"])^2 *~~(x * 10^3)))


pdf(file.path(figures_path, "abias-amse.pdf"), width = 12, height = 6)
print(plot_bias + plot_mse + plot_layout(axes = "collect"))
dev.off()

