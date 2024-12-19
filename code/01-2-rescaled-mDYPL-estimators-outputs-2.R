supp_path = "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")
figure <- "2"

library("dplyr")
library("ggplot2")
library("colorspace")
library("ggpp")
library("patchwork")

source(file.path(supp_path, "code/methods/plot-with-insets.R"))
load(file.path(results_path, paste0('rescaled-mDYPL-estimates-figure-', figure, '.rda')))

alpha_fac <- 500
p_size <- 0.5
## exp_h <- 5 * 200 * 2
vp_h <- 0.14
## exp_w <- exp_h * sqrt(2)
## exp_ratio <- exp_h / exp_w
vp_w <- vp_h * 1.1
cols <- qualitative_hcl(4)

p_ests <- plot_with_insets(ests, pt, type = "estimate_and_truth", N = n, col = cols[c(1, 3, 4)]) +
    labs(y = expression(gamma), title = "MDYPL estimator")

p_rescaled_ests <- plot_with_insets(rescaled_ests, pt, type = "estimate_and_truth", N = n, col = cols[c(1, 3, 4)]) +
    labs(y = expression(gamma), title = "rescaled MDYPL estimator")

pdf(file.path(figures_path, paste0("mdypl-vs-truth-", figure,".pdf")), width = 9, height = 4)
print(p_ests + p_rescaled_ests)
dev.off()
