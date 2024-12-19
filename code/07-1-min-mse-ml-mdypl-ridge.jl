n_cores = 10

using Distributed
using RCall
using CSV

addprocs(n_cores)

@everywhere begin
    using ProgressMeter
    using SharedArrays
    using Optim, JLD2, LineSearches, NonlinearSolve, Distributions, NLsolve
end

@everywhere begin
    supp_path = "."
    results_path = joinpath(supp_path, "results")
    include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
    include(joinpath(home_dir, "code", "methods", "AMP_Ridge.jl"))
end

@everywhere using .mDYPL
@everywhere using .AMP_Ridge

@rput supp_path results_path
R"""
library("parallel")
source(file.path(supp_path, "code/methods/compute-pt.R"))
""";

kappas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
gammas = [1., 2.5, 5, 7.5, 10, 12.5, 15.0]
kg = reduce(vcat, Base.product(kappas, gammas))
kg = reduce(vcat, transpose(map(v -> collect(v), kg)))
nkg = size(kg)[1]

## Compute phase transition
@rput kg n_cores
R"""
library("parallel")
kg <- as.data.frame(kg)
colnames(kg) <- c("kappa", "gamma")
set.seed(123)
ns <- 200000
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
kg$kappa_diff <- mclapply(seq.int(nrow(kg)), function(j) {
    gamma <- kg$gamma[j]
    kappa <- kg$kappa[j]
    compute_pt(beta0 = 0, gamma, ncores = 1, XZU =  xzu)$kappa - kappa
}, mc.cores = ncores) |> unlist()
mle_exists <- kg$kappa_diff > 0
"""
@rget mle_exists

grid_size = 100
alphas = vcat(range(0.01, 0.9, 2 * grid_size),
              range(0.901, 0.99, 3 * grid_size))
lambdas = vcat(range(4.5, 0.1, 2 * grid_size),
               range(0.101, 0.001, 3 * grid_size))

len = length(alphas)
report = len

## ML
ML_pars = SharedArray{Float64}(nkg, 3)
ML_res = SharedArray{Float64}(nkg, 3)
ML_pars[:, :, :] .= NaN
ML_res[:, :, :] .= NaN
s_ML = [1., 2., 1.] ./ 5
@showprogress @distributed for i in 1:nkg
    if mle_exists[i]
        κ = kg[i, 1]
        γ = kg[i, 2]
        o = Vector{Float64}(undef, 3)
        start_ML = s_ML
        ML = solve_mDYPL_SE(κ, 1., γ,
                            verbose = false,
                            method = NewtonRaphson(),
                            x_init = start_ML,
                            abstol = 1e-10,
                            reltol = 1e-10)
        mDYPL.mDYPL_SE!(o, ML, κ, γ, 1.,  use_eta = false)
        if (any(abs.(o) .> 0.0001))
            start_ML = [0.5, 0.5, 1.5]
            ML = solve_mDYPL_SE(κ, 1., γ, 
                                verbose = false,
                                method = TrustRegion(),
                                x_init = start_ML,
                                abstol = 1e-10,
                                reltol = 1e-10)
            mDYPL.mDYPL_SE!(o, ML, κ, γ, 1.,  use_eta = false)
        end
        ML_pars[i, :] = ML
        ML_res[i, :] = o
        println("\nML | ",
                i, " / ", nkg, " | ",
                "κ = ", round(κ, digits = 2), ", γ = ", round(γ, digits = 2), " | ",
                round(maximum(abs.(ML_res[i, :])), digits = 12))
    else
        ML_pars[i, :] .= NaN
        ML_res[i, :] .= NaN
    end    
end

maximum(filter(!isnan, abs.(ML_res))) 


## mDYPL
MDYPL_pars = SharedArray{Float64}(nkg, len, 3)
MDYPL_res = SharedArray{Float64}(nkg, len, 3)
MDYPL_pars[:, :, :] .= NaN
MDYPL_res[:, :, :] .= NaN
s_MDYPL = [1., 2., 1.] ./ 2
@showprogress @distributed for i in 1:nkg
    κ = kg[i, 1]
    γ = kg[i, 2]
    o = Vector{Float64}(undef, 3)
    for j in 1:len
        if j == 1
            start_MDYPL = s_MDYPL
        else
            start_MDYPL = MDYPL_pars[i, j - 1, :]
        end
        α = alphas[j]
        MDYPL = solve_mDYPL_SE(κ, α, γ, 
                               verbose = false,
                               method = TrustRegion(),
                               x_init = start_MDYPL,
                               abstol = 1e-10,
                               reltol = 1e-10)
        mDYPL.mDYPL_SE!(o, MDYPL, κ, γ, α, use_eta = false)
        if (any(abs.(o) .> 0.0001))
            start_MDYPL = [0.5, 0.5, 1.5]
            MDYPL = solve_mDYPL_SE(κ, α, γ,
                                   verbose = false,
                                   method = TrustRegion(),
                                   x_init = start_MDYPL,
                                   abstol = 1e-10,
                                   reltol = 1e-10)
            mDYPL.mDYPL_SE!(o, MDYPL, κ, γ, α, use_eta = false)
        end
        MDYPL_pars[i, j, :] = MDYPL
        MDYPL_res[i, j, :] = o
        if j % report == 0
            println("\nmDYPL | ",
                    i, " / ", nkg, " | ",
                    j, " / ", len, " | ",             
                    "κ = ", round(κ, digits = 2), ", γ = ", round(γ, digits = 2),
                    ", α = ", round(α, digits = 4), " | ",
                    round(maximum(abs.(MDYPL_res[i, 1:j, :])), digits = 12))
        end
    end
end

maximum(filter(!isnan, abs.(MDYPL_res))) 


## ridge
ridge_pars = SharedArray{Float64}(nkg, len, 3)
ridge_res = SharedArray{Float64}(nkg, len, 3)
ridge_pars[:, :, :] .= NaN
ridge_res[:, :, :] .= NaN
s_ridge = [1., 2., 1.] ./ 5
@showprogress @distributed for i in 1:nkg
    κ = kg[i, 1]
    γ = kg[i, 2]
    o = Vector{Float64}(undef, 3)
    for j in 1:len
        if j == 1
            start_ridge = s_ridge
        else
            start_ridge = ridge_pars[i, j - 1, :]
        end
        λ = lambdas[j]
        ridge = find_params_ridge_nonlinearsolve(κ, γ, λ, 
                                                 verbose = false,
                                                 method = NewtonRaphson(),
                                                 x_init = start_ridge,
                                                 abstol = 1e-10,
                                                 reltol = 1e-10)
        AMP_Ridge.eq_bin!(o, ridge, κ, γ, λ)
        if (any(abs.(o) .> 0.0001))
            start_ridge = [0.5, 0.5, 1.5]
            ridge = find_params_ridge_nonlinearsolve(κ, γ, λ, 
                                                     verbose = false,
                                                     method = TrustRegion(),
                                                     x_init = start_ridge,
                                                     abstol = 1e-10,
                                                     reltol = 1e-10)
            AMP_Ridge.eq_bin!(o, ridge, κ, γ, λ)
        end
        ridge_pars[i, j, :] = ridge
        ridge_res[i, j, :] = o
        if j % report == 0
            println("\nmDYPL | ",
                    i, " / ", nkg, " | ",
                    j, " / ", len, " | ",             
                    "κ = ", round(κ, digits = 2), ", γ = ", round(γ, digits = 2),
                    ", λ = ", round(λ, digits = 4), " | ",
                    round(maximum(abs.(ridge_res[i, 1:j, :])), digits = 12))
        end
    end
end

maximum(filter(!isnan, abs.(ridge_res)))

kgambs = Matrix{Float64}(undef, nkg * len, 7)
for j in 1:nkg
    kgambs[(j - 1) * len .+ (1:len), :] = hcat(reduce(vcat, map(i -> kg[j, :]', 1:len)),
                                               alphas,
                                               MDYPL_pars[j, :, :],
                                               maximum(abs.(MDYPL_res[j, :, :]), dims = 2))
end


kglmbs = Matrix{Float64}(undef, nkg * len, 7)
for j in 1:nkg
    kglmbs[(j - 1) * len .+ (1:len), :] = hcat(reduce(vcat, map(i -> kg[j, :]', 1:len)),
                                               lambdas,
                                               ridge_pars[j, :, :],
                                               maximum(abs.(ridge_res[j, :, :]), dims = 2))
end

kgmbs = hcat(kg, ML_pars, maximum(abs.(ML_res), dims = 2))

@rput kgmbs
R"""
library("ggplot2")
library("dplyr")
kgmbs <- data.frame(kgmbs)
names(kgmbs) <- c("kappa", "gamma", "mu", "b", "sigma", "linf")
"""

@rput kgambs
R"""
kgambs <- data.frame(kgambs)
names(kgambs) <- c("kappa", "gamma", "alpha", "mu", "b", "sigma", "linf")
"""

@rput kglmbs
R"""
kglmbs <- data.frame(kglmbs)
names(kglmbs) <- c("kappa", "gamma", "lambda", "mu", "b", "sigma", "linf")
"""


out_mDYPL = joinpath(results_path, "state_evolution_mDYPL.rda")
out_ridge = joinpath(results_path, "state_evolution_ridge.rda")
out_ML = joinpath(results_path, "state_evolution_ML.rda")

R"""
save(kgambs, file = $out_mDYPL)
save(kglmbs, file = $out_ridge)
save(kgmbs, file = $out_ML)
"""











