n_cores = 10

using Distributed
using RCall

addprocs(n_cores)

@everywhere begin
    using Random, Optim, NonlinearSolve, InvertedIndices
    using SharedArrays, ProgressMeter
    supp_path = "."
    results_path = joinpath(supp_path, "results")
end

@everywhere begin
    include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
end

@everywhere using .mDYPL

@rput supp_path results_path
R"""
library("parallel")
source(file.path(supp_path, "code/methods/compute-pt.R"))
""";

kappas = [0.1, 0.5]
alphas = [1, 3/4, 1/2, 1/4]
gamma = sqrt(5)

kga = vec(collect(Base.product(kappas, gamma, alphas)))
kga = reduce(vcat, map(x -> hcat(x[1], x[2], x[3]), kga))
kga = vcat(kga, [0.1, gamma, 1 / (1 + 0.1)]', [0.5, gamma, 1 / (1 + 0.5)]')


## check when mle exists
@rput gamma n_cores
R"""
ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
pt_point <- compute_pt(gamma_grid = gamma, ncores = n_cores, XZU = xzu)
"""
@rget pt_point

mle_exists = kga[:, 1] .< pt_point[:, "kappa"]
n_settings = size(kga)[1]

@everywhere begin
    function scale_beta(beta, n, gamma)
        return sqrt(n) * gamma * beta / sqrt(beta' * beta)
    end
end

R = 1000
n = 2000

I_5 = 1:5
I_50 = 1:50

loglik_full = SharedArray{Float64}(n_settings, R)
loglik_nest_5 = SharedArray{Float64}(n_settings, R)
loglik_nest_50 = SharedArray{Float64}(n_settings, R)

converged_full = SharedArray{Bool}(n_settings, R)
converged_nest_5 = SharedArray{Bool}(n_settings, R)
converged_nest_50 = SharedArray{Bool}(n_settings, R)

## Get R * n_settings unique seeds
Random.seed!(123)
seeds = Matrix{Int}(undef, n_settings, R)
while length(unique(seeds)) != R * n_settings
    seeds[:, :] = floor.(Int64, rand(n_settings, R) .* 10000000)
end

@showprogress @distributed for (r, i) in collect(Iterators.product(1:R, 1:n_settings))
    kappa = kga[i, 1]
    a = kga[i, 3]
    gamma = kga[i, 2]
    p = Int(ceil(n * kappa))
    Random.seed!(seeds[i, r])
    X = randn(n, p) ./ sqrt(n)
    beta = vcat(zeros(Int(p / 2)), ones(Int(p / 2)))
    beta = scale_beta(beta, n, gamma)
    mu = 1.0 ./ (1.0 .+ exp.(.- X * beta))
    y = rand(n) .< mu
    fit_full = logistic_mDYPL(y, X, a)
    beta_DY_coefs = Optim.minimizer(fit_full)
    fit_nest_5 = logistic_mDYPL(y, X[:, Not(I_5)], a)
    fit_nest_50 = logistic_mDYPL(y, X[:, Not(I_50)], a)
    converged_full[i, r] = Optim.converged(fit_full)
    converged_nest_5[i, r] = Optim.converged(fit_nest_5)
    converged_nest_50[i, r] = Optim.converged(fit_nest_50)
    loglik_full[i, r] = Optim.minimum(fit_full)
    loglik_nest_5[i, r] = Optim.minimum(fit_nest_5)
    loglik_nest_50[i, r] = Optim.minimum(fit_nest_50)
end

plr_5 = -2 * (loglik_full .- loglik_nest_5);
plr_50 = -2 * (loglik_full .- loglik_nest_50);

consts = Matrix{Float64}(undef, n_settings, 3)
for i in 1:n_settings
    kappa = kga[i, 1]
    gamma = kga[i, 2]
    a = kga[i, 3]
    if a == 1.0 && !mle_exists[i]
        consts[i, :] .= NaN
    else
        consts[i, :] = solve_mDYPL_SE(kappa, a, gamma, use_eta = false,
                                      verbose = true,
                                      x_init = [0.5, 0.4, 1.0],
                                      method = TrustRegion())       
    end
end



@rput plr_5 plr_50 kga consts
R"""
colnames(kga) <- c("kappa", "gamma", "alpha")
kga <- as.data.frame(kga)
kga <- kga |>
   transform(alpha_fac = as.character(MASS::fractions(alpha)),
             kappa_fac = paste("kappa ==", kappa))
kga[kga$alpha_fac %in% c("10/11", "2/3"), "alpha_fac"] <- "1 / (1 + kappa)"
kga$alpha_fac <- factor(kga$alpha_fac, levels = c("1", "3/4", "1/2", "1/4", "1 / (1 + kappa)"), ordered = TRUE)
levels(kga$alpha_fac) <- paste("alpha ==", levels(kga$alpha_fac))
reshape_mat <- function(mat, method) {
    data.frame(value = c(mat), method = method)
}
adj_plr_5 <- plr_5 / (kga[, 1] * consts[, 3]^2 / consts[, 2])
lr_5 <- rbind(reshape_mat(plr_5, method = "PLR"),
              reshape_mat(adj_plr_5, method = "rescaled PLR"))
lr_5 <- data.frame(lr_5, kga, df = 5)
adj_plr_50 <- plr_50 / (kga[, 1] * consts[, 3]^2 / consts[, 2])
lr_50 <- rbind(reshape_mat(plr_50, method = "PLR"),
              reshape_mat(adj_plr_50, method = "rescaled PLR"))
lr_50 <- data.frame(lr_50, kga, df = 50)
"""

R"""
save.image(file.path(results_path, 'qqplots-rescaled-plr.rda'))
"""
