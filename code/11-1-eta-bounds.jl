using Distributions, LinearAlgebra, RCall, Optim, NonlinearSolve, ProgressMeter, Random, SharedArrays 
supp_path = "." 
results_path = joinpath(supp_path, "results")
figures_path = joinpath(supp_path, "figures")
include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))

using .mDYPL

function get_oracle_constants(kappa, gamma, alpha, init)
    result = fill(NaN, 3)
    consts = solve_mDYPL_SE(kappa, alpha, gamma; use_eta = false,
                            numeric = false,
                            method = TrustRegion(),
                            abstol = 1e-9, reltol = 1e-9,
                            x_init = init)
    mDYPL.mDYPL_SE!(result, consts, kappa, gamma, alpha; use_eta = false)
    Dict("constants" => consts,
         "results" => result)
end

function get_upsilon_gamma(kappa, init; alpha_grid, gamma_grid = 0.0001:0.1:20)
    n_gamma = length(gamma_grid)
    upsilon = Vector{Float64}(undef, n_gamma)
    constants = Matrix{Float64}(undef, n_gamma, 3)
    results = Matrix{Float64}(undef, n_gamma, 3)
    for j in 1:n_gamma
        consts = get_oracle_constants(kappa, gamma_grid[j], alpha_grid[j], init)
        init = consts["constants"]
        constants[j, :] = consts["constants"]
        results[j, :] = consts["results"]
        upsilon[j] = sqrt(init[1]^2 * gamma_grid[j]^2 + kappa * init[3]^2)
        # println(j, " κ = ", kappa, " γ = ", gamma_grid[j], " ", round(maximum(abs.(consts["results"])), digits = 10))
    end
    ## alpha, kappa, gamma, upsilon, mu, b, sigma, f_mu, f_b, f_sigma
    hcat(alpha_grid, ones(n_gamma) .* kappa, gamma_grid, upsilon, constants, results)
end


se_init = [1.2, 1.2, 3.]

gamma_grid = 0.0001:0.1:15
kappas = 0.1:0.1:0.9
alphas = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]

settings = vec(collect(Base.product(kappas, alphas)))
settings = mapreduce(x -> hcat(x[1], x[2]), vcat, settings)
n_settings = size(settings)[1]

ug = Vector{Matrix}(undef, n_settings)
@showprogress Threads.@threads for k in 1:n_settings
    kappa = settings[k, 1]
    alpha = settings[k, 2]
    alpha_grid = alpha * ones(length(gamma_grid))
    ug[k] = get_upsilon_gamma(kappa, se_init, gamma_grid = gamma_grid,
                              alpha_grid = alpha_grid)
end

ug_adapt = Vector{Matrix}(undef, length(kappas))
alpha_adapt = @. exp(gamma_grid / 2) / (1 + exp(gamma_grid / 2))
@showprogress Threads.@threads for k in 1:length(kappas)
    ug_adapt[k] = get_upsilon_gamma(kappas[k], se_init, gamma_grid = gamma_grid,
                                    alpha_grid = alpha_adapt)
end


@rput ug ug_adapt results_path
R"""
library("ggplot2")
ug <- do.call("rbind", ug)
ug <- data.frame(ug)
ug_adapt <- do.call("rbind", ug_adapt)
ug_adapt <- data.frame(ug_adapt)
names(ug_adapt) <- names(ug) <- c("alpha", "kappa", "gamma", "upsilon",
                                  "mu", "b", "sigma", "f_mu", "f_b", "f_sigma")
save(ug, ug_adapt, file = file.path(results_path, "upsilon-gamma.rda"))
"""



