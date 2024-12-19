n_cores = 10

using Distributed

addprocs(n_cores)

@everywhere begin 
    using Distributions, LinearAlgebra, RCall, Optim, NonlinearSolve, ProgressMeter, Random, SharedArrays
    supp_path = "." 
    results_path = joinpath(supp_path, "results")
    include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
end

@everywhere begin
    using .mDYPL
end

@everywhere begin 
    function simu_X(n, kappa)
        p = floor(Int, kappa * n) 
        Sigma = fill(0.5, p, p) + 0.5 * I(p)
        mvnorm = MvNormal(zeros(p), Sigma)
        Matrix(rand(mvnorm, n)')
    end

    function get_beta(n, kappa, gamma)
        p = floor(Int, kappa * n) 
        Sigma = fill(0.5, p, p) + 0.5 * I(p)
        inds = 1:floor(Int, 5)
        not_inds = setdiff(1:p, inds)
        beta = fill(0., p)
        beta[inds] .= 1.
            mult = beta' * Sigma * beta
        beta / sqrt(mult) * gamma        
    end

    function simu_y(X, beta)
        n = size(X)[1]
        eta = X * beta
        mu = @. 1.0 / (1.0 + exp(- eta))
        rand(n) .< mu
    end

    function eta_lasso(X, y)
        @rput X y
        R"""
      library("glmnet")
      cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)
      best_lambda <- cv_fit$lambda.min
      best_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)
      coefs <- coef(best_model) 
      lasso_eta <-  predict(best_model, newx = X, type = "link")
        """
        @rget lasso_eta
    end

    function cls(X, y, eta)
        ym = 2 * y .- 1.
        XTX = transpose(X) * X
        invXTX = inv(XTX)
        Px = X * invXTX * transpose(X)
        betaHat0 = invXTX * transpose(X) * ym
        eta_norm = dot(eta, eta)
        eta_tanh = tanh.(eta ./2)
        cHat = transpose(eta) * eta_tanh / eta_norm
        Phat = Px - eta * transpose(eta) / eta_norm
        deltaHat = invXTX * transpose(X) * Phat * eta_tanh        
        (betaHat0 - deltaHat) / cHat
    end

    function mdypl(X, y, alpha)
        Optim.minimizer(logistic_mDYPL(y, X, alpha))
    end

    function get_oracle_constants(kappa, gamma, alpha, init)
        result = fill(NaN, 3)
        consts = solve_mDYPL_SE(kappa, alpha, gamma; use_eta = false,
                                numeric = false,
                                abstol = 1e-9, reltol = 1e-9,
                                x_init = init)
        mDYPL.mDYPL_SE!(result, consts, kappa, gamma, alpha; use_eta = false)
        Dict("constants" => consts,
             "results" => result)
    end

    function get_constants(X, y, alpha, beta_init, init)
        result = fill(NaN, 3)
        upsilon_hat = get_eta(y, X, alpha; beta_init = beta_init)
        d = size(X)
        k = d[2] / d[1]
        consts = solve_mDYPL_SE(k, alpha, upsilon_hat; use_eta = true,
                                numeric = false,
                                method = TrustRegion(),
                                verbose = false,
                                abstol = 1e-9, reltol = 1e-9,
                                x_init = init)
        mDYPL.mDYPL_SE!(result, consts, k, upsilon_hat, alpha; use_eta = true)
        Dict("constants" => consts,
             "results" => result,
             "upsilon" => upsilon_hat)
    end

    function estimate_tau(X)
        n, p = size(X) 
        chol = cholesky(X' * X) 
        function jrss(j)
            e_j = zeros(p)
            e_j[j] = 1
            v_j = chol \ e_j
            RSS_j = 1 / v_j[j]
            sqrt(RSS_j  / (n - p + 1))
        end
        map(jrss, 1:p)
    end  
    
    function LB_z_stat(estimate, beta0, X, eta)
        sigma_hat = dot(eta, tanh.(eta ./ 2)) / dot(eta, eta)
        n = length(eta)
        buff = inv(X'X) * X' * (sqrt.(1. .- tanh.(eta ./2).^2) .* I(n))
        sigma_hat * (estimate .- beta0) ./ sqrt.(diag(buff * buff'))
    end
   
end

R = 1000
gamma = 3.
alpha = 0.95
ns = [400, 800, 1200, 1600, 2000]
kappas = [0.2, 0.5]

settings = vec(collect(Base.product(ns, kappas)))
settings = mapreduce(x -> hcat(x[1], x[2]), vcat, settings)
n_settings = size(settings)[1]

Random.seed!(123)
seeds = Matrix{Int}(undef, n_settings, R)
while length(unique(seeds)) != R * n_settings
    seeds[:, :] = floor.(Int64, rand(n_settings, R) .* 1000000000)
end

se_init = [1.2, 1.2, 3.]

for i in 1:n_settings 

    nobs = Int64(settings[i, 1])
    kappa = settings[i, 2]
    
    Random.seed!(1)
    X = simu_X(nobs, kappa)
    beta_true = get_beta(nobs, kappa, gamma)
    p = length(beta_true)
    eta_true = Vector(X * beta_true)
    se_oracle = get_oracle_constants(kappa, gamma, alpha, se_init)

    y = SharedArray{Float64}(nobs, R)
    eta_hat = SharedArray{Float64}(nobs, R)
    @showprogress @distributed for r in 1:R
        Random.seed!(seeds[i, r])
        y[:, r] = simu_y(X, beta_true)
        eta_hat[:, r] = eta_lasso(X, y[:, r])
    end

    simu_inds = 1:R
    beta_mdypl = mapreduce(r -> mdypl(X, y[:, r], alpha), hcat, simu_inds)

    ##
    ## 3 consts, 3 SE functions, 1 upsilon
    consts = SharedArray{Float64}(7, R)
    @showprogress @distributed for r in 1:R
        se = get_constants(X, y[:, r], alpha, beta_mdypl[:, r], se_oracle["constants"])
        consts[1:3, r] = se["constants"]
        consts[4:6, r] = se["results"]
        consts[7, r] = se["upsilon"]
    end

    se_converged = mapslices(r -> maximum(abs.(r)), consts[4:6, :], dims = 1) .< 1e-4   
    sum(se_converged)
    
    ## MDYPL
    tau = sqrt((p + 1) / (2 * p))
    beta_mdypl_rescaled_oracle = beta_mdypl ./ se_oracle["constants"][1]
    rescaled_z_oracle = (beta_mdypl .- se_oracle["constants"][1] * zeros(p)) / (se_oracle["constants"][3] / (sqrt(nobs) * tau))
    beta_mdypl_rescaled = beta_mdypl ./ consts[1, :]'
    tau_hat = estimate_tau(X)
    rescaled_z = (sqrt(nobs) .* tau_hat) .* beta_mdypl ./ (consts[3, :])'

    ## CLS
    beta_cls_oracle = mapreduce(r -> cls(X, y[:, r], eta_true), hcat, simu_inds)
    beta_cls = mapreduce(r -> cls(X, y[:, r], eta_hat[:, r]), hcat, simu_inds)
    z_cls_oracle = mapreduce(r -> LB_z_stat(beta_cls_oracle[:, r], zeros(p), X, eta_true), hcat, simu_inds)
    z_cls = mapreduce(r -> LB_z_stat(beta_cls[:, r], zeros(p), X, eta_hat[:, r]), hcat, simu_inds)

    @rput R nobs kappa gamma results_path
    @rput z_cls z_cls_oracle rescaled_z rescaled_z_oracle se_converged consts
    @rput beta_mdypl beta_cls beta_cls_oracle beta_mdypl_rescaled beta_mdypl_rescaled_oracle beta_true
R"""
library(dplyr)
p <- length(beta_true)
methds <- c("CLS", "CLS [O]")
n_methds <- length(methds)
stats <-  c("CLS", "CLS [O]")
n_stats <- length(stats)
estimates <- data.frame(estimate = c(beta_cls, beta_cls_oracle),
                        method = rep(methds, each = R * p),
                        truth = beta_true,
                        parameter = rep(1:p, R * n_methds),
                        sample = rep(rep(1:R, each = p), n_methds),
                        kappa = kappa, gamma = gamma, N = nobs, p = p,
                        converged = TRUE)
z_statistics <- data.frame(z = c(z_cls, z_cls_oracle),
                           statistic = rep(stats, each = R * p),
                           truth = beta_true,
                           parameter = rep(1:p, R * n_stats),
                           sample = rep(rep(1:R, each = p), n_stats),
                           kappa = kappa, gamma = gamma, N = nobs, p = p,
                           converged = TRUE)
save(estimates, 
     file = file.path(results_path, paste0("cLS-MDYPL-comparison-estimates-cLS-n=", nobs, "-kappa=", kappa, ".rda")))
save(z_statistics,
     file = file.path(results_path, paste0("cLS-MDYPL-comparison-statistics-cLS-n=", nobs, "-kappa=", kappa, ".rda")))
""";

R"""
p <- length(beta_true)
methds <- c("rescaled MDYPL", "rescaled MDYPL [O]")
n_methds <- length(methds)
stats <-  c("rescaled MDYPL", "rescaled MDYPL [O]")
n_stats <- length(stats)
estimates <- data.frame(estimate = c(beta_mdypl_rescaled, beta_mdypl_rescaled_oracle),
                        method = rep(methds, each = R * p),
                        converged = TRUE,
                        truth = beta_true,
                        parameter = rep(1:p, R * n_methds),
                        sample = rep(rep(1:R, each = p), n_methds),
                        kappa = kappa, gamma = gamma, N = nobs, p = p)
estimates <- estimates |> mutate(converged = ifelse(method == "rescaled MDYPL", rep(se_converged, each = unique(p)), TRUE))
z_statistics <- data.frame(z = c(rescaled_z, rescaled_z_oracle),
                           statistic = rep(stats, each = R * p),
                           truth = beta_true,
                           parameter = rep(1:p, R * n_stats),
                           sample = rep(rep(1:R, each = p), n_stats),
                           kappa = kappa, gamma = gamma, N = nobs, p = p)
z_statistics <- z_statistics |> mutate(converged = ifelse(statistic == "rescaled MDYPL", rep(se_converged, each = unique(p)), TRUE))
constants <- data.frame(t(consts))
names(constants) <- c("mu", "b", "sigma", "f_mu", "f_b", "f_sigma", "upsilon")
constants <- constants |> mutate(sample = 1:R,
                                 kappa = kappa, gamma = gamma, N = nobs, p = p)
save(estimates, constants,
     file = file.path(results_path, paste0("cLS-MDYPL-comparison-estimates-MDYPL-n=", nobs, "-kappa=", kappa, ".rda")))
save(z_statistics, 
     file = file.path(results_path, paste0("cLS-MDYPL-comparison-statistics-MDYPL-n=", nobs, "-kappa=", kappa, ".rda")))
""";
    
end


## Save oracle constats
se_oracle = SharedArray{Float64}(n_settings, 3)
@showprogress @distributed for i in 1:n_settings 
    nobs = Int64(settings[i, 1])
    kappa = settings[i, 2]
    X = simu_X(nobs, kappa)
    beta_true = get_beta(nobs, kappa, gamma)
    p = length(beta_true)
    eta_true = Vector(X * beta_true)
    se_oracle[i, :] = get_oracle_constants(kappa, gamma, alpha, se_init)["constants"]
end

@rput se_oracle settings results_path
R"""
se_oracle <- data.frame(se_oracle)
names(se_oracle) <- c("mu", "b", "sigma")
settings <- data.frame(settings)
names(settings) <- c("n", "kappa")
settings <- cbind(settings, se_oracle)
settings <- settings |> transform(root_aMSE = sqrt(2 / (n + 1 / kappa)) * sigma / mu)
save(settings, file = file.path(results_path, "cLS-MDYPL-comparison-settings.rda"))
"""

# ## Indices of zero coefficients
# ind_zero = findall(beta_true .== 0)
## Indices of non-zero coefficients
# ind_nonzero = findall(beta_true .!= 0)
# what = ind_zero[1]
# @rput rescaled_z_oracle z_cls_oracle rescaled_z z_cls what failed
# R"""
# par(mfrow = c(2, 2))
# hist(2 * pnorm(-abs(rescaled_z_oracle[what, ])), prob = TRUE, main = "rescaled z [oracle]", xlab = "z")
# abline(h = 1)
# hist(2 * pnorm(-abs(z_cls_oracle[what, ])), prob = TRUE, main = "cls z [oracle]", xlab = "z")
# abline(h = 1)
# hist(2 * pnorm(-abs(rescaled_z[what, ])), prob = TRUE, main = "rescaled z", xlab = "z")
# abline(h = 1)
# hist(2 * pnorm(-abs(z_cls[what, ])), prob = TRUE, main = "cls z", xlab = "z")
# abline(h = 1)
# c(mean(abs(rescaled_z_oracle[what, ]) > qnorm(0.975)),
#       mean(abs(rescaled_z[what, ]) > qnorm(0.975)),
#       mean(abs(z_cls_oracle[what, ]) > qnorm(0.975)),
#       mean(abs(z_cls[what,]) > qnorm(0.975)))
# """

# root_aMSE(beta_cls, beta_true)
# root_aMSE(beta_cls_oracle, beta_true)
# root_aMSE(beta_mdypl, beta_true)
# root_aMSE(beta_mdypl_rescaled, beta_true)
# root_aMSE(beta_mdypl_rescaled_oracle, beta_true)
# sqrt(2 * (se_oracle["constants"][3]^2 / se_oracle["constants"][1]^2) / (nobs + 1 / kappa))

# bias(beta_cls, beta_true)["bias"][[ind_nonzero[1], ind_zero[1]]]
# bias(beta_cls_oracle, beta_true)["bias"][[ind_nonzero[1], ind_zero[1]]]
# bias(beta_mdypl, beta_true)["bias"][[ind_nonzero[1], ind_zero[1]]]
# bias(beta_mdypl_rescaled, beta_true)["bias"][[ind_nonzero[1], ind_zero[1]]]
# bias(beta_mdypl_rescaled_oracle, beta_true)["bias"][[ind_nonzero[1], ind_zero[1]]]
         
