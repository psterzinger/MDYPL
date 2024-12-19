## Useful functions
get_variables <- function(terms, covariate_names) {
    vars <- NULL
    for (te in terms) {
       vars <- c(vars, grep(te, covariate_names, value = TRUE))
    }
    vars
}

simu_fun <- function(coefs, mm) {
    probs <- drop(plogis(mm %*% coefs))
    list(y = as.numeric(runif(nrow(mm)) < probs), X = mm)
}

adjust_response <- function(y, alpha) {
    alpha * y + (1 - alpha) / 2
}

get_loglik <- function(obj) {
    with(obj, sum(y * log(fitted.values)) + sum((1 - y) * log(1 - fitted.values)))
}

solve_SE_intercept <- function(kappa, alpha, eta, beta0hat, start, warm = FALSE,
                               code_path = NULL) {
    aux <- data.frame(kappa = kappa, alpha = alpha, eta = eta, beta0hat = beta0hat)
    file_name <- paste0(c(sample(letters, 15, replace = TRUE),
                          sample(1:10, 5, replace = TRUE)), collapse = "")
    file_name <- paste0(file_name, ".rda", collapse = "")
    save(aux, start, warm, file = file_name)
    system(paste("julia", file.path(code_path, "compute_constants_intercept.jl"), file_name))
    load(file_name)
    file.remove(file_name)
    list(consts = consts, funs = funs)
}


out_name <- function(terms, nsimu, beta0, gamma) {
    base_name <- paste0(paste0(terms, collapse = "+"), "-simu")
    paste0(base_name,
           "-R=", nsimu,
           "-beta0=", beta0,
           "-gammmasq=", gamma^2,
           ".rda")
}
