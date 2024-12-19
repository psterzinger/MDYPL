supp_path <- "."
results_path = file.path(supp_path, "results")
figures_path = file.path(supp_path, "figures")
code_path <- file.path(supp_path, "code")
data_path <- file.path(supp_path, "data")

library("lmtest")
library("detectseparation")
library("ggplot2")
library("colorspace")
library("patchwork")
library("glmnet")

source(file.path(code_path, "methods", "mfeat-functions.R"))

chosen_digit <- 7

full_terms <- c("fou", "kar")
nest_terms <- c("fou")

## Constract a data frame with response and all features
feat_types <- c("fac", "fou", "kar", "mor", "pix", "zer")
mfeat <- as.list(numeric(length(feat_types)))
names(mfeat) <- feat_types
for (feat_type in feat_types) {
    cdat <- read.table(file.path(data_path, paste0("mfeat-", feat_type)))
    colnames(cdat) <- seq.int(ncol(cdat))
    mfeat[[feat_type]] <- cdat
}
feats <- sapply(mfeat, ncol)
feats <- data.frame(type = names(feats), number = feats)
mfeat_data <- data.frame(digit = factor(rep(0:9, each = 200)), do.call("cbind", mfeat))

chosen_digit <- as.character(chosen_digit)

## Add response
mfeat_data <- mfeat_data |>
    transform(y = digit == chosen_digit)

all_covariates <- names(mfeat_data)[!(names(mfeat_data) %in% c("digit", "y"))]
all_formula <- formula(paste("y ~ ", paste(all_covariates, collapse = " + ")))
Xall <- model.matrix(all_formula, data = mfeat_data)

## Check for any aliased parameters and remove them from the model matrix
qrXall <- qr(Xall)
if (qrXall$rank < ncol(Xall)) {
    aliased <- qrXall$pivot[seq.int(qrXall$rank + 1, ncol(Xall))]
    Xall <- Xall[, -aliased]
}
all_covariates <- colnames(Xall)
nobs_all <- nrow(Xall)

## Training and test data
set.seed(123) ## same as simulation
train_id <- sample(1:nobs_all, nobs_all * 0.5, replace = FALSE)
X_train <- Xall[train_id, ]
X_test <- Xall[-train_id, ]
y_train <- mfeat_data[train_id, "y"]
y_test <- mfeat_data[-train_id, "y"]
nobs <- nrow(X_train)

## Center X's
X_train[, -1] <- scale(X_train[, -1], center = TRUE, scale = FALSE)
X_test[, -1] <- scale(X_test[, -1], center = TRUE, scale = FALSE)

## Get model matrices for ful and nested model
full_X <- X_train[, full_vars <- get_variables(full_terms, all_covariates), drop = FALSE]
nest_X <- X_train[, nest_vars <- get_variables(nest_terms, all_covariates), drop = FALSE]
full_p <- ncol(full_X)
nest_p <- ncol(nest_X)

## Characteristics
kappa <- ncol(full_X) / nrow(full_X)
alpha <- 1 / (1 + kappa)
degrees <- full_p - nest_p

## Adjusted responses
ya_train <- adjust_response(y_train, alpha)
full_train_data <- data.frame(y = ya_train, full_X)
nest_train_data <- data.frame(y = ya_train, nest_X)

full_fit <- glm(y ~ ., data = full_train_data, family = binomial())
nest_fit <- glm(y ~ ., data = nest_train_data, family = binomial())
full_ll <- get_loglik(full_fit)
nest_ll <- get_loglik(nest_fit)
plr = 2 * (full_ll - nest_ll)

## Get SLOE estimate of eta
mu <- fitted(full_fit)
v <- mu * (1 - mu)
h <- hatvalues(full_fit)
S <- full_fit$linear.predictors - (ya_train - mu) / v * h / (1 - h)
eta <- sd(S)

results <- data.frame(plr = plr,
                      eta = eta,
                      n = nrow(full_X),
                      p = full_p,
                      kappa = kappa,
                      alpha = alpha,
                      df = degrees,
                      beta0hat = coef(full_fit)[1])

starting_values <- expand.grid(mu = c(0.5),
                               b = c(1., 0.5),
                               sigma = c(3., 1.),
                               beta0 = c(-2, 0))
starting_values <- as.matrix(t(starting_values))
SE_sol <- solve_SE_intercept(results$kappa, results$alpha, results$eta, results$beta0hat,
                             starting_values, code_path = file.path(code_path, "methods"))
results <- data.frame(results, SE_sol$consts, SE_sol$funs) |>
    transform(rescaled_plr = b * plr / (kappa * sigma^2),
              gammahat = sqrt(eta^2 - kappa * sigma^2) / mu)




## Outputs

## Fits on original data
full_train_data_o <- data.frame(y = y_train, full_X)
nest_train_data_o <- data.frame(y = y_train, nest_X)
full_fit_o <- update(full_fit, data = full_train_data_o)
nest_fit_o <- update(nest_fit, data = full_train_data_o)
full_fit_sep <- update(full_fit_o, method = detect_separation)
nest_fit_sep <- update(nest_fit_o, method = detect_separation)
full_vs_nest_o <- lrtest(nest_fit_o, full_fit_o)

## statistics, p-values, infinite estimates
tab <- data.frame(test = c("LR", "PLR", "rescaled PLR"),
                  statistic = c(full_vs_nest_o[2, "Chisq"],
                                results$plr,
                                results$rescaled_plr),
                  pvalue = c(full_vs_nest_o[2, "Pr(>Chisq)"],
                             pchisq(results$plr, degrees, lower.tail = FALSE),
                             pchisq(results$rescaled_plr, degrees, lower.tail = FALSE)),
                  full_coefs_inf = c(full_fit_sep$outcome, FALSE, FALSE),
                  nest_coefs_inf = c(nest_fit_sep$outcome, FALSE, FALSE))
tab

## Fitted values and predictions on test set
dsize = 2

train_sample <- full_train_data |>
    transform(digit = mfeat_data[train_id, "digit"],
              y = y_train)
train_sample <- train_sample[order(train_sample$digit), ]

y_test <- mfeat_data[-train_id, "y"]
digits_test <- mfeat_data[-train_id, "digit"]
test_sample <- cbind(data.frame(y = y_test, digit = digits_test), X_test)
test_sample <- test_sample[order(test_sample$digit), ]


train_sample1 <- train_sample |>
    transform(eta = predict(full_fit, newdata = train_sample, type = "link")) |>
    transform(probs = plogis(eta),
              estimator = "mDYPL",
              instance = 1:nrow(train_sample),
              data_subset = "training")
train_sample2 <- train_sample1 |>
    transform(eta = results$intercept + (eta - coef(full_fit)[1]) / results$mu) |>
    transform(probs = plogis(eta),
             estimator = "rescaled mDYPL",
             instance = 1:nrow(train_sample),
             data_subset = "training")

test_sample1 <- test_sample |>
    transform(eta = predict(full_fit, newdata = test_sample, type = "link")) |>
    transform(probs = plogis(eta),
              estimator = "mDYPL",
              instance = 1:nrow(test_sample),
              data_subset = "test")
test_sample2 <- test_sample1 |>
    transform(eta = results$intercept + (eta - coef(full_fit)[1]) / results$mu) |>
    transform(probs = plogis(eta),
              estimator = "rescaled mDYPL",
              instance = 1:nrow(test_sample),
              data_subset = "test")

keep <- c("y", "digit", "probs", "estimator", "instance", "data_subset")
all_data <- rbind(rbind(train_sample1, train_sample2)[keep],
                  rbind(test_sample1, test_sample2)[keep]) |>
    transform(data_subset = factor(data_subset, levels = c("training", "test"), ordered = TRUE))


### LASSO
if (FALSE) {
    elnet <- cv.glmnet(full_X, y_train, family = "binomial")
    oo <- order(mfeat_data[-train_id, "digit"])
    fit_test <- predict(elnet, newx = X_test[oo, full_vars], type = "response")
    oo <- order(mfeat_data[train_id, "digit"])
    fit_train <- predict(elnet, newx = X_train[oo, full_vars], type = "response")
    test_lasso <- test_sample2 |> transform(probs = fit_test,
                                            estimator = "lasso",
                                            instance = 1:nrow(test_sample),
                                            data_subset = "test")
    train_lasso <- train_sample2 |> transform(probs = fit_train,
                                              estimator = "lasso",
                                              instance = 1:nrow(train_sample),
                                              data_subset = "training")
    all_data <- rbind(all_data, test_lasso[keep], train_lasso[keep]) |>
        transform(data_subset = factor(data_subset, levels = c("training", "test"), ordered = TRUE))
}

probs_plot <- ggplot(all_data) +
    geom_hline(aes(yintercept = 0.5), lty = 2, col = "grey") +
    geom_text(aes(x = instance, y = probs, label = digit, col = digit != chosen_digit),
              size = dsize) +
    facet_grid(data_subset ~ estimator) +
    theme_minimal() +
    labs(x = "Observation", y = "Probability") +
    scale_color_grey(guide = "none") +
    labs(title = "fou + kar") +
    lims(y = c(0, 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())



## nested model - mDYPL
train_sample3 <- train_sample |>
    transform(eta = predict(nest_fit, newdata = train_sample, type = "link")) |>
    transform(probs = plogis(eta),
              estimator = "mDYPL",
              instance = 1:nrow(train_sample),
              data_subset = "training")

test_sample3 <- test_sample |>
    transform(eta = predict(nest_fit, newdata = test_sample, type = "link")) |>
    transform(probs = plogis(eta),
              estimator = "mDYPL",
              instance = 1:nrow(test_sample),
              data_subset = "test")

all_data_nest <- rbind(train_sample3[keep],
                       test_sample3[keep]) |>
    transform(data_subset = factor(data_subset,
                                   levels = c("training", "test"),
                                   ordered = TRUE))

probs_plot_nest <- ggplot(all_data_nest) +
    geom_hline(aes(yintercept = 0.5), lty = 2, col = "grey") +
    geom_text(aes(x = instance, y = probs, label = digit, col = digit != chosen_digit),
              size = dsize) +
    facet_grid(data_subset ~ estimator) +
    theme_minimal() +
    labs(x = "Observation", y = "Probability") +
    scale_color_grey(guide = "none") +
    labs(title = "fou") +
    lims(y = c(0, 1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y = element_blank())


pdf(file.path(figures_path, "case-study-probs.pdf"), width = 6, height = 4)
print(probs_plot_nest + probs_plot + plot_layout(widths = c(1, 3), axes = "collect"))
dev.off()


