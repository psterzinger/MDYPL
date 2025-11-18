chosen_digit <- 7
n_cores <- 10
nsimu <- 500
supp_path <- "."

results_path = file.path(supp_path, "results")
code_path <- file.path(supp_path, "code")
data_path <- file.path(supp_path, "data")

## Define full and nested models
full_terms <- c("fou", "kar")
nest_terms <- c("fou")

### Combinations of gamma and beta0 to attempt
gb0 <- expand.grid(gamma = sqrt(c(1, 2, 4, 8, 16)),
                   beta0 = c(-3, -2, -1, 0))

library("parallel")
source(file.path(code_path, "methods", "mfeat-functions.R"))

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

## modelling
chosen_digit <- as.character(chosen_digit)

## Add response
mfeat_data <- mfeat_data |>
    transform(y = digit == chosen_digit)

all_covariates <- names(mfeat_data)[!(names(mfeat_data) %in% c("digit", "y"))]
all_formula <- formula(paste("y ~ ", paste(all_covariates, collapse = " + ")))
Xall <- model.matrix(all_formula, data = mfeat_data)

## Check for aliased parameters and remove them from the model matrix
qrXall <- qr(Xall)
if (qrXall$rank < ncol(Xall)) {
    aliased <- qrXall$pivot[seq.int(qrXall$rank + 1, ncol(Xall))]
    Xall <- Xall[, -aliased]
}
all_covariates <- colnames(Xall)
nobs_all <- nrow(Xall)

## Train data
set.seed(123) ## same as case study
train_id <- sample(1:nobs_all, 0.5 * nobs_all, replace = FALSE)
X <- Xall[train_id, ]
nobs <- nrow(X)

X[, -1] <- scale(X[, -1], center = TRUE, scale = FALSE)

## Simulation study

full_X <- X[, full_vars <- get_variables(full_terms, all_covariates), drop = FALSE]
nest_X <- X[, nest_vars <- get_variables(nest_terms, all_covariates), drop = FALSE]
full_p <- ncol(full_X)
nest_p <- ncol(nest_X)

## Characteristics
kappa <- ncol(full_X) / nrow(full_X)
alpha <- 1 / (1 + kappa)
degrees <- full_p - nest_p


for (j in seq.int(nrow(gb0))) {
    gamma <- gb0[j, "gamma"]
    beta0 <- gb0[j, "beta0"]
    cat(j, "/", nrow(gb0), ": Computing for gamma^2 =", gamma^2, "beta0 =", beta0, "\n")

    ## Truth
    set.seed(333)
    full_coefs <- structure(rep(0, full_p), names = full_vars)
    full_coefs[nest_vars] <- rnorm(nest_p)
    full_coefs <- c(beta0, gamma * full_coefs / sd(full_X %*% full_coefs))

    ## Simulate data sets
    simu_data <- replicate(nsimu, simu_fun(full_coefs, cbind(1, full_X)), simplify = FALSE)
    fresults0 <- mclapply(seq.int(nsimu), function(k) {
        cat(k, "/", nsimu, "\n")
        y_o <- simu_data[[k]]$y
        y_adj <- adjust_response(y_o, alpha)
        cX <- simu_data[[k]]$X[, -1]
        full_train_data <- data.frame(y_adj = y_adj, cX)
        nest_train_data <- data.frame(y_adj = y_adj, cX[, nest_vars])
        ## Compute likelihood ratio statistic
        full_md <- glm(y_adj ~ ., family = binomial(), data = full_train_data)
        nest_md <- glm(y_adj ~ ., family = binomial(), data = nest_train_data)
        full_ll <- get_loglik(full_md)
        nest_ll <- get_loglik(nest_md)
        plr = 2 * (full_ll - nest_ll)
        ## Get SLOE estimate of eta
        mu <- fitted(full_md)
        v <- mu * (1 - mu)
        h <- hatvalues(full_md)
        S <- full_md$linear.predictors - (y_adj - mu) / v * h / (1 - h)
        eta <- sd(S)
        list(coefficients = data.frame(estimate = coef(full_md),
                                       truth = full_coefs,
                                       parameter = names(coef(full_md)),
                                       sample = k),
             plr = data.frame(plr = plr,
                              eta = eta,
                              beta0hat = coef(full_md)[1],
                              sample = k))
    }, mc.cores = n_cores)
    results <- do.call("rbind", lapply(results0, "[[", "plr"))
    coefs <- do.call("rbind", lapply(results0, "[[", "coefficients"))

    results$df <- degrees
    results$n <- nobs
    results$p <- full_p
    results$alpha <- alpha
    results$kappa <- kappa

    ## Get constants
    starting_values <- expand.grid(mu = c(0.3, 0.5),
                                   b = c(2., 1., 0.5),
                                   sigma = c(2, 2, 3.),
                                   beta0 = c(-3, -2, 0, 2))
    starting_values <- as.matrix(t(starting_values))
    SE_sol <- solve_SE_intercept(results$kappa, results$alpha, results$eta, results$beta0hat,
                                 starting_values, code_path = file.path(code_path, "methods"))

    results <- data.frame(results, SE_sol$consts, SE_sol$funs) |>
        transform(rescaled_plr = b * plr / (kappa * sigma^2),
                  gammahat = sqrt(eta^2 - kappa * sigma^2) / mu,
                  beta0 = beta0,
                  gamma = gamma,
                  kappa = kappa) |>
        transform(converged = (abs(mu_fun) < 1e-08) &
                      (abs(b_fun) < 1e-08) &
                      (abs(sigma_fun) < 1e-08) &
                      (abs(intercept_fun) < 1e-08))

    keep <- c("estimate", "truth", "rescaled_estimate", "parameter",
              "converged", "gamma", "beta0", "sample")
    coefs <- coefs |>
        merge(results, by = "sample") |>
        transform(rescaled_estimate = ifelse(parameter == "(Intercept)",
                                             intercept,
                                             estimate / mu)) |>
        subset(select = keep)

    save(results, coefs,
         file = file.path(results_path, out_name(full_terms, nsimu, beta0, gamma)))

}
