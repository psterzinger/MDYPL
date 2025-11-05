generate_unique_seeds <- function(n, maxint = n^2) {
    out <- round(runif(n, 0, maxint))
    test <- any(dinds <- duplicated(out))
    while (test) {
        out[dinds] <- round(runif(sum(dinds), 0, .Machine$integer.max/2))
        test <- any(dinds <- duplicated(out))
    }
    out
}
