supp_path <- "."
figures_path <- file.path(supp_path, "figures")

library("brglm2")
library("detectseparation")
library("ggplot2")
library("patchwork")

data("MultipleFeatures", package = "brglm2")

## Center the fou.* and kar.* features in training and test sets
vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
train_id <- which(MultipleFeatures$training)
MultipleFeatures[train_id, vars] <- scale(MultipleFeatures[train_id, vars], scale = FALSE)
MultipleFeatures[-train_id, vars] <- scale(MultipleFeatures[-train_id, vars], scale = FALSE)

## Compute the MDYPL fits
kappa <- length(vars) / sum(MultipleFeatures$training)
full_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))
nest_vars <- grep("fou", vars, value = TRUE)
nest_fm <- formula(paste("I(digit == 7) ~", paste(nest_vars, collapse = " + ")))

full_mdypl <- glm(full_fm, data = MultipleFeatures, family = binomial(), subset = training,
                  method = mdyplFit, alpha = 1 / (1 + kappa))
nest_mdypl <- update(full_mdypl, nest_fm)
full_ml <- glm(full_fm, data = MultipleFeatures, family = binomial(), subset = training)
nest_ml <- update(full_ml, nest_fm)

## Detect separation
full_sep <- update(full_ml, method = detect_separation)
nest_sep <- update(nest_ml, method = detect_separation)

full_vs_nest_lr <- anova(nest_ml, full_ml)
full_vs_nest_plr <- plrtest(nest_mdypl, full_mdypl)
full_vs_nest_rplr <- plrtest(nest_mdypl, full_mdypl, hd_correction = TRUE)

## Statistics
tab0 <- data.frame(test = c("LR", "PLR", "rescaled PLR"),
                  statistic = c(full_vs_nest_lr[2, "Deviance"],
                                full_vs_nest_plr[2, "Deviance"],
                                full_vs_nest_rplr[2, "Deviance"]),
                  pvalue = c(full_vs_nest_lr[2, "Pr(>Chi)"],
                             full_vs_nest_plr[2, "Pr(>Chi)"],
                             full_vs_nest_rplr[2, "Pr(>Chi)"]),
                  full_infinite = c(full_sep$outcome, FALSE, FALSE),
                  nest_infinite = c(nest_sep$outcome, FALSE, FALSE))
tab0 |> transform(statistic = round(statistic, 2),
                  pvalue = round(pvalue, 2))

## Plots
## MDYPL
MultipleFeatures$prob_full <- predict(full_mdypl, newdata = MultipleFeatures,
                                      type = "response")
MultipleFeatures$prob_nest <- predict(nest_mdypl, newdata = MultipleFeatures,
                                      type = "response")
MultipleFeatures$estimator <- "MDYPL"
MultipleFeatures <- MultipleFeatures |>
    transform(type = ifelse(training, "training", "test"))
## Rescaled MDYPL
MultipleFeaturesR <- MultipleFeatures
MultipleFeaturesR$estimator <- "rescaled MDYPL"
full_summ <- summary(full_mdypl, hd_correction = TRUE)
X_full <- model.matrix(full_fm, data = MultipleFeaturesR)
MultipleFeaturesR$prob_full <- plogis(drop(X_full %*% coef(full_summ)[, "Estimate"]))
MultipleFeaturesR <- MultipleFeaturesR |>
    transform(type = ifelse(training, "training", "test"))

dsize <- 2
probs_nest_mdypl <- ggplot(MultipleFeatures) +
    geom_hline(aes(yintercept = 0.5), lty = 2, col = "grey") +
    geom_text(aes(x = 1:nrow(MultipleFeatures), y = prob_nest, label = digit, col = digit != 7),
              size = dsize) +
    facet_grid(type ~ estimator) +
    theme_minimal() +
    labs(x = "Observation", y = "Probability") +
    scale_color_grey(guide = "none") +
    labs(title = "fou") +
    lims(y = c(0, 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y = element_blank())

probs_full_mdypl <- ggplot(MultipleFeatures) +
    geom_hline(aes(yintercept = 0.5), lty = 2, col = "grey") +
    geom_text(aes(x = 1:nrow(MultipleFeatures), y = prob_full, label = digit, col = digit != 7),
              size = dsize) +
    facet_grid(type ~ estimator) +
    theme_minimal() +
    labs(x = "Observation", y = "Probability") +
    scale_color_grey(guide = "none") +
    labs(title = "fou + kar") +
    lims(y = c(0, 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y = element_blank())

probs_full_r_mdypl <- ggplot(MultipleFeaturesR) +
    geom_hline(aes(yintercept = 0.5), lty = 2, col = "grey") +
    geom_text(aes(x = 1:nrow(MultipleFeaturesR), y = prob_full, label = digit, col = digit != 7),
              size = dsize) +
    facet_grid(type ~ estimator) +
    theme_minimal() +
    labs(x = "Observation", y = "Probability") +
    scale_color_grey(guide = "none") +
    labs(title = "fou + kar") +
    lims(y = c(0, 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

pdf(file.path(figures_path, "case-study-probs.pdf"), width = 6, height = 4)
print(probs_nest_mdypl + probs_full_mdypl + probs_full_r_mdypl +  plot_layout(axes = "collect"))
dev.off()
