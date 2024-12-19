supp_path = "."
results_path = joinpath(supp_path, "results")
n_cores = 10
figure = "2"

using Random, Optim, NonlinearSolve, RCall
include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
using .mDYPL

@rput supp_path results_path
R"""
library("parallel")
source(file.path(supp_path, "code/methods/compute-pt.R"))
""";

function fill_beta(p,b)
    beta = Vector{Float64}(undef,p) 
    pp = length(b) 
    for i in eachindex(beta) 
        beta[i] = b[(i-1)%pp+1]
    end 
    sort(beta)
end 

function scale_beta(beta, n, gamma)
    return sqrt(n) * gamma * beta / sqrt(beta' * beta)
end

kappas = (0.1, 0.2, 0.4, 0.6, 0.8, 0.9) 
gammas = (2, 6, 10, 14, 18)
full_grid = [(0.0,0.0) for i in  Iterators.product(gammas,kappas)]
for i in eachindex(gammas)
    gamma = gammas[i]
    for j in eachindex(kappas)
        kappa = kappas[j]
        full_grid[i,j] = (kappa,gamma) 
    end 
end 
grid_inds = vcat(1, 11, 21, 7, 17, 3, 13, 23, 9, 19, 5, 15, 25)
kg_pairs = vcat((0.2, sqrt(0.9)), (sqrt(.9) / 22.4, 0.2 * 22.4), full_grid[grid_inds][2:end])

if figure == "1"
    n = 2000
    b = [-3, -3/2, 0, 3/2, 3]
elseif figure == "2"
    n = 1000
    b = vcat([-10, 10], zeros(6))
end

reps = 10
betas = [Float64[] for i in eachindex(kg_pairs)] 
beta_true = [Float64[] for i in  eachindex(kg_pairs)] 
params = [Float64[] for i in  eachindex(kg_pairs)]
counter = 1
for (kappa, gamma) in kg_pairs 
    println("κ:$kappa, ", "γ:$gamma")
    p = floor(Int64, n * kappa)
    beta = fill_beta(p, b)
    beta = scale_beta(beta, n, gamma)
    beta_true[counter] = beta
    a = 1 / (1 + kappa)
    for k in 1:reps 
        Random.seed!((counter - 1) * length(kg_pairs) + k) 
        X = randn(n, p) / sqrt(n) 
        mu = 1.0 ./ (1.0 .+ exp.(.- X*beta)) 
        y = rand(n) .< mu 
        beta_DY = Optim.minimizer(logistic_mDYPL(y, X, a; beta_init = beta))
        if k == 1
            betas[counter] = beta_DY
        else 
            betas[counter] .= betas[counter] .+ beta_DY 
        end
    end
    params[counter] = solve_mDYPL_SE(kappa, a, gamma, use_eta = false, method = TrustRegion())
    counter += 1
end 
betas = betas ./ reps
rescaled_betas = betas ./ [x[1] for x in params]
kg_pairs = reduce(vcat, map(x -> hcat(x[1], x[2]), kg_pairs))


@rput betas rescaled_betas kg_pairs n_cores beta_true n figure
R"""
ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
ga <- c(0.001, 0.01, seq(0, 20, length = 200))
pt <- compute_pt(gamma_grid = ga, ncores = n_cores, XZU = xzu)
""";

R"""
nbeta = lengths(betas)
ests <- data.frame(
   estimate = unlist(betas),
   kappa = rep(kg_pairs[, 1], nbeta),
   gamma = rep(kg_pairs[, 2], nbeta),
   truth = unlist(beta_true),
   parameter = unlist(sapply(nbeta, function(i) 1:i)),
   mle_exists = TRUE,
   method = "mDYPL"
)
rescaled_ests <- within(ests, {
   estimate = unlist(rescaled_betas)
})
""";

R"""
save(n, pt, ests, rescaled_ests, 
     file = file.path(results_path, paste0('rescaled-mDYPL-estimates-figure-', figure, '.rda')))
"""


