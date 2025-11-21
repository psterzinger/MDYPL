fit_mdypl <- function(X, y, alpha = 1, start = rep(0, ncol(X)),
                      eps_f = 1e-14, eps_g = 1e-14, maxit = 1000) {
    ya <- alpha * y + (1 - alpha) / 2
    res <- fastLR(X, ya, start, eps_f, eps_g, maxit = maxit)
    res$pl <- sum(ya * log(fitted(res)) + (1 - ya) * log(1 - fitted(res)))
    res
}
