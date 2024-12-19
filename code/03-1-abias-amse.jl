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

kappas = .0075:.005:.999
gammas = .025:0.1:20
n_kappas = length(kappas)
n_gammas = length(gammas)

consts = SharedArray{Float64}(n_kappas, n_gammas, 3)
result = SharedArray{Float64}(n_kappas, n_gammas, 3)
start = [0.5, 0.5, 1.]
@showprogress @distributed for ig in 1:n_gammas
    gamma = gammas[ig]
    for ik in 1:n_kappas
        kappa = kappas[ik]
        o = Vector{Float64}(undef, 3)
        if ik == 1
            init = start
        else
            init = consts[ik - 1, ig, :]
        end
        consts[ik, ig, :] = solve_mDYPL_SE(kappa, 1 / (1 + kappa), gamma,
                                           use_eta = false,
                                           verbose = false,
                                           x_init = init,
                                           abstol = 1e-10,
                                           reltol = 1e-10,
                                           method = TrustRegion())
        mDYPL.mDYPL_SE!(o, consts[ik, ig, :], kappa, gamma, 1 / (1 + kappa),
                        use_eta = false)
        result[ik, ig, :] = o
        println("MDYPL | ",
                "i_gamma: ", ig, " / ", n_gammas, " | ",
                "i_kappa: ", ik, " / ", n_kappas, " | ",
                "κ = ", round(kappa, digits = 2), ", γ = ", round(gamma, digits = 2), " | ",
                round(maximum(abs.(result[1:ik, ig, :])), digits = 12))
    end
end

@rput gammas n_cores
R"""
ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
pts <- compute_pt(gamma_grid = gammas, ncores = n_cores, XZU = xzu)
"""

@rput kappas gammas consts
R"""
df <- NULL
dd <- dim(consts)
n_kappas <- dd[1]
n_gammas <- dd[2]
for (i in seq.int(n_kappas)) {
   cdat <- consts[i,,]
   cdat <- as.data.frame(cdat)
   names(cdat) <- c("mu", "b", "sigma")
   cdat$kappa <- kappas[i]
   cdat$gamma <- gammas
   df <- rbind(df, cdat)
}
"""

R"""
save(df, pts, file = file.path(results_path, 'abias-amse.rda'))
"""
