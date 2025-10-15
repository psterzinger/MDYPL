supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")

load(file.path(results_path, "upsilon-gamma.rda"))

library("ggplot2")

pdf(file.path(figures_path, "upsilon-gamma.pdf"), width = 8, height = 2)
ggplot(ug) +
   geom_line(aes(gamma, upsilon, color = alpha, group = alpha)) +
   geom_line(data = ug_adapt, aes(gamma, upsilon), col = "grey", linewidth = 1) +
   facet_grid(~ kappa, labeller = label_bquote(cols = kappa == .(kappa))) +
   labs(x = expression(gamma), y = expression(upsilon)) +
   scale_color_continuous(type = "viridis") +
   coord_cartesian(y = range(ug$upsilon)) +
   theme_minimal()
dev.off()

files <- dir(results_path, "cLS-MDYPL-comparison-estimates-MDYPL-n", full.names = TRUE)
consts <- NULL
for (f in files) {
    load(f)
    consts <- rbind(consts, constants)
}

bounds <- ug |>
    subset(kappa %in% c(0.2, 0.5) & alpha == 0.95) |>
    aggregate(upsilon ~ kappa + alpha, FUN = \(x) c(min = min(x), max = max(x))) |>
    merge(unique(consts[c("kappa", "N")]), by = "kappa")
bounds$minimum <- bounds$upsilon[, "min"]
bounds$maximum <- bounds$upsilon[, "max"]


pdf(file.path(figures_path, "cLS-MDYPL-upsilon.pdf"), width = 8, height = 3)
ggplot(consts) +
    geom_histogram(aes(upsilon), bins = 50) +
    geom_rect(data = bounds, aes(xmin = minimum, xmax = maximum, ymin = 0, ymax = Inf), alpha = 0.2) +
    facet_grid(kappa ~ N, labeller = label_bquote(cols = N == .(N), rows = kappa == .(kappa))) +
    labs(x = expression(upsilon)) +
    theme_minimal()
dev.off()
