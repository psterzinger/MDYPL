## Phase transition curves
compute_pt <- function(beta0 = 0,
                       gamma_grid = seq(0, 20, length = 100),
                       nsimu = 1000000,
                       ncores = 10,
                       XZU = NULL) {
    ## Theorem 2.1 of DOI: 10.1214/18-AOS1789
    if (is.null(XZU)) {
        X <- rnorm(nsimu)
        U <- runif(nsimu)
        Z <- rnorm(nsimu)
    } else {
        X <- XZU$X
        U <- XZU$U
        Z <- XZU$Z
    }
    kappa <- function(gamma) {
        gamma0 <- sqrt(gamma^2 - beta0^2)
        Y <- -1 + 2 * (plogis(beta0 + gamma0 * X) > U)
        V <- X * Y
        obj <- function(ts) {
            mean(pmax(ts[1] * Y + ts[2] * V - Z, 0)^2)
        }
        optim(c(0, 0), obj, method = "BFGS")$value
    }
    kappa <- future_sapply(1:length(gamma_grid), function(k) {
        kappa(gamma = gamma_grid[k])
    })
    data.frame(kappa = kappa, gamma = gamma_grid)
}
