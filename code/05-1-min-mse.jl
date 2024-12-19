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

kappas = .05:.0125:.975
n_kappas = length(kappas)
gammas = .05:0.25:10.0
n_gammas = length(gammas)

kg = reduce(vcat, Base.product(kappas, gammas))
kg = reduce(vcat, transpose(map(v -> collect(v), kg)))
n_kg = size(kg)[1]

grid_size = 20
alphas = vcat(range(0.01, 0.95, 3 * grid_size),
              range(0.95, 0.99, 2 * grid_size))
n_alphas = length(alphas)

consts = SharedArray{Float64}(n_kg, n_alphas, 3);
result = SharedArray{Float64}(n_kg, n_alphas, 3);
start = [0.5, 1.0, 1.0];

start = [0.9, 0.1, 0.1];

@showprogress @distributed for i in 1:n_kg
    kappa = kg[i, 1]
    gamma = kg[i, 2]
    o = Vector{Float64}(undef, 3)
    for ia in 1:n_alphas
        if ia == 1
            init = start
         else
            init = consts[i, ia - 1, :] .* [1.0, 1.1, 1.04]
        end
        alpha = alphas[ia]
        consts[i, ia, :] = solve_mDYPL_SE(kappa, alpha, gamma,
                                          use_eta = false,
                                          verbose = false,
                                          x_init = init,
                                          abstol = 1e-10,
                                          reltol = 1e-10,
                                          method = TrustRegion())
        mDYPL.mDYPL_SE!(o, consts[i, ia, :], kappa, gamma, alpha,
                        use_eta = false)
        if (any(abs.(o) .> 0.00001))
            consts[i, ia, :] = solve_mDYPL_SE(kappa, alpha, gamma,
                                              use_eta = false,
                                              verbose = false,
                                              x_init = init,
                                              abstol = 1e-10,
                                              reltol = 1e-10,
                                              method = NewtonRaphson(),
                                              numeric = true)
        end
        if (any(abs.(o) .> 0.00001))
            consts[i, ia, :] = solve_mDYPL_SE(kappa, alpha, gamma,
                                              use_eta = false,
                                              verbose = false,
                                              x_init = start,
                                              abstol = 1e-10,
                                              reltol = 1e-10,
                                              method = TrustRegion())
            mDYPL.mDYPL_SE!(o, consts[i, ia, :], kappa, gamma, alpha,  use_eta = false)
        end       
        result[i, ia, :] = o
    end   
    println("\nMDYPL | ",
            "i = ", i, " / ", n_kg, " | ",
            "κ = ", round(kappa, digits = 3),
            ", γ = ", round(gamma, digits = 3), " | ",
            round(maximum(abs.(result[i, :, :])), digits = 12))
end

## Another attempt to non-converged
not_converged = findall(abs.(result) .> 1e-06)
if length(not_converged) > 0 
    inds = unique(reduce(hcat, map(x -> [x[1], x[2]], not_converged)), dims = 2)
    all(inds[2, :] .== 1.)

    for k in 1:size(inds)[2]
        i = inds[1, k]
        ia = inds[2, k]
        kappa = kg[i, 1]
        gamma = kg[i, 2]
        o = Vector{Float64}(undef, 3)
        init = consts[i, ia + 1, :] ./ [1.0, 1.1, 1.04]
        alpha = alphas[ia]
        consts[i, ia, :] = solve_mDYPL_SE(kappa, alpha, gamma,
                                          use_eta = false,
                                          verbose = true,
                                          x_init = init,
                                          abstol = 1e-10,
                                          reltol = 1e-10,
                                          method = TrustRegion())
        mDYPL.mDYPL_SE!(o, consts[i, ia, :], kappa, gamma, alpha,
                        use_eta = false)
        result[i, ia, :] = o
        println("\nMDYPL | ",
                "i = ", i, " / ", n_kg, " | ",
                "κ = ", round(kappa, digits = 3),
                ", γ = ", round(gamma, digits = 3), " | ",
                round(maximum(abs.(result[i, :, :])), digits = 12))
    end
end

## Another attempt to iny mus
not_converged = findall(abs.(consts[:, :, 1]) .< 1e-02)
if length(not_converged) > 0 
    inds = unique(reduce(hcat, map(x -> [x[1], x[2]], not_converged)), dims = 2)
    all(inds[2, :] .== 1.)

    for k in 1:size(inds)[2]
        i = inds[1, k]
        ia = inds[2, k]
        kappa = kg[i, 1]
        gamma = kg[i, 2]
        o = Vector{Float64}(undef, 3)
        init = consts[i, ia + 1, :] ./ [1.0, 1.1, 1.04]
        alpha = alphas[ia]
        consts[i, ia, :] = solve_mDYPL_SE(kappa, alpha, gamma,
                                          use_eta = false,
                                          verbose = true,
                                          x_init = init,
                                          abstol = 1e-10,
                                          reltol = 1e-10,
                                          method = TrustRegion())
        mDYPL.mDYPL_SE!(o, consts[i, ia, :], kappa, gamma, alpha,
                        use_eta = false)
        result[i, ia, :] = o
        println("\nMDYPL | ",
                "i = ", i, " / ", n_kg, " | ",
                "κ = ", round(kappa, digits = 3),
                ", γ = ", round(gamma, digits = 3), " | ",
                round(maximum(abs.(result[i, :, :])), digits = 12))
    end
end


alphas_min_mse = map(i -> alphas[argmin(consts[i, :, 3].^2 ./ consts[i, :, 1].^2)], 1:n_kg)

## check when mle exists
@rput gammas n_cores
R"""
ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
pts <- compute_pt(gamma_grid = gammas, ncores = n_cores, XZU = xzu)
"""

@rput alphas_min_mse kg consts result
R"""
kg <- data.frame(kg)
names(kg) <- c("kappa", "gamma")
df <- data.frame(alpha = alphas_min_mse, kg)
"""

R"""
save(df, pts, consts, result, file = file.path(results_path, 'alpha-min-mse.rda'))
"""
