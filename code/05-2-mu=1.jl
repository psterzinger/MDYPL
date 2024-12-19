n_cores = 10

using Distributed
using RCall

using Random, Optim, NonlinearSolve, InvertedIndices
using SharedArrays, ProgressMeter
supp_path = "."
results_path = joinpath(supp_path, "results")

include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))

using .mDYPL

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

consts = SharedArray{Float64}(n_kg, 3);
result = SharedArray{Float64}(n_kg, 3);
start = [0.5, 1.0, 1.0];

@showprogress for i in 1:n_kg
    kappa = kg[i, 1]
    gamma = kg[i, 2]
    o = Vector{Float64}(undef, 3)
    if i == 1
        init = start
    elseif (i - 1) % n_kappas == 0
        id = Int(n_kappas * ((i - 1) / n_kappas - 1) + 1)
        init = consts[id, :]
    else
        init = consts[i - 1, :] .* [1.0, 1.2, 1.02]
    end
    ## consts has alpha, b, sigma
    consts[i, :] = solve_mDYPL_SE(kappa, NaN, gamma; 
                                  use_eta = false, 
                                  numeric = true,
                                  method = TrustRegion(), 
                                  verbose = false,
                                  mu_fixed = true,
                                  x_init = init,
                                  mu_target = 1.0,
                                  abstol = 1e-10, 
                                  reltol = 1e-10) 
    mDYPL.mDYPL_SE!(o, [1.0, consts[i, 2], consts[i, 3]] ,
                    kappa, gamma, consts[i, 1],
                    use_eta = false)
    result[i, :] = o   
    println("\nMDYPL | ",
            "i = ", i, " / ", n_kg, " | ",
            "κ = ", round(kappa, digits = 3),
            ", α = ", round(consts[i, 1], digits = 4),
            ", γ = ", round(gamma, digits = 3), " | \n",
            "solu: ", round.(consts[i, :], digits = 2), " | \n",
            "init: ", round.(init, digits = 2), " | ",
            round(maximum(abs.(result[i, :])), digits = 12))
end

## check when mle exists
@rput gammas n_cores
R"""
ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
pts <- compute_pt(gamma_grid = gammas, ncores = n_cores, XZU = xzu)
"""

@rput consts kg pts result
R"""
df <- data.frame(consts, kg)
names(df) <- c("alpha", "b", "sigma", "kappa", "gamma")
"""

R"""
save(df, pts, file = file.path(results_path, 'mu=1.rda'))
"""
