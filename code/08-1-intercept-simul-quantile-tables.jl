n_cores = 10

using Distributed, DataFrames, PrettyTables, CSV

addprocs(n_cores) 

@everywhere begin 
    using Random, Optim, LinearAlgebra, NonlinearSolve, Statistics, SharedArrays, Distributions, ProgressMeter
    supp_path = "." 
    results_path = joinpath(supp_path, "results")

    include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
    @eval using .mDYPL
    intercepts = 0.5:.5:2.5
    kappas = .2:.2:.8  
    gamma = sqrt(5) 
    reps = 5000
    n = 2000
    dfs = 10 
end 

pars_true = fill(NaN, length(kappas), length(intercepts), 4)
for (i, kappa) in enumerate(kappas) 
    for (j, intercept) in enumerate(intercepts) 
        if i == j == 1 
            x_init = vcat(solve_mDYPL_SE(kappa, 1 / (1 + kappa), gamma; use_eta = false, method = TrustRegion(), verbose = false, abstol = 1e-9, reltol = 1e-9), .65)
        elseif j == 1 
            x_init = pars_true[i - 1, j, :] 
        else
            x_init = pars_true[i, j - 1, :] 
        end 
        pars_true[i, j, :] = solve_mDYPL_SE(kappa, 1 / (1 + kappa), gamma, intercept; use_eta = false, method = TrustRegion(), verbose = true, abstol = 1e-9, reltol = 1e-9, x_init = x_init, iterations = 5000)  
    end 
end 
@everywhere pars_true = $pars_true

# Simul 
null_T_stats = SharedArray{Float64}(length(kappas), length(intercepts), reps) 
T_stats = SharedArray{Float64}(length(kappas), length(intercepts), reps) 
llr = SharedArray{Float64}(length(kappas), length(intercepts), reps)
for i in eachindex(kappas) 
    kappa = kappas[i]
    for j in eachindex(intercepts) 
        intercept = intercepts[j] 
        p = floor(Int, kappa * n) 
        Sigma = I(p) / n 

        Random.seed!((i-1) * length(intercepts) + j) 
        inds = rand(1:p, floor(Int, p / 4))   
        not_inds = setdiff(1:p,inds)
        # ind = rand(inds) + 1
        # ind2 = rand(not_inds) + 1 
        ind = inds[1] + 1
        ind2 = not_inds[1] + 1
        # llr_inds = sample(not_inds, dfs, replace = false) .+ 1 
        llr_inds = not_inds[1:dfs] .+ 1 

        beta = fill(0.,p) 
        beta[inds] .= 1.
        mult = beta' * Sigma * beta
        beta = beta / sqrt(mult) * gamma 
        beta = vcat(intercept, beta) 

        mvnorm = MvNormal(zeros(p), Sigma) 
       
        tau = 1 / sqrt(n) 
        null_tau = 1 / sqrt(n) 
        X = fill(1., n, p + 1)
        X_nested = fill(1., n, p + 1 - dfs) 
        y = fill(NaN,n) 
        res = fill(NaN,4)
        eta = X * beta
        mu = 1.0 ./ (1.0 .+ exp.(.-eta))
        beta_est = fill(NaN, p + 1)
        a = 1 / (1 + kappa) 
        @everywhere begin 
            beta = $beta 
            eta = $eta 
            mu = $mu 
            X = $X 
            kappa = $kappa 
            inds = $inds 
            not_inds = $not_inds
            y = $y 
            beta_est = $beta_est 
            kappa = $kappa 
            n = $n 
            mvnorm = $mvnorm
            ind = $ind 
            ind2 = $ind2 
            tau = $tau
            null_tau = $null_tau
            llr_inds = $llr_inds
            X_nested = $X_nested
            a = $a
        end
        println("Setting $(((i-1) * length(intercepts) + j)) / $(length(kappas) * length(intercepts)): κ = $kappa, β₀ = $intercept")
        @showprogress @distributed for k in 1:reps
            Random.seed!((i-1) * length(intercepts) * reps + (j - 1) * reps + k)
            @views X[:, 2:end] = rand(mvnorm, n)' 
            @views X_nested .= X[:, setdiff(1:(p + 1), llr_inds)] 
            eta .= X * beta
            mu .= 1.0 ./ (1.0 .+ exp.(.-eta))
            y .= rand(n) .< mu
            
            mDYPL_full = logistic_mDYPL(y, X, a) 
            mDYPL_nested = logistic_mDYPL(y, X_nested, a)
            beta_est .= Optim.minimizer(mDYPL_full)

            llr[i,j,k] = 2 * (Optim.minimum(mDYPL_nested) - Optim.minimum(mDYPL_full)) * pars_true[i, j, 2] / (pars_true[i, j, 3]^2 * kappa) 
            null_T_stats[i,j,k] = sqrt(n) * beta_est[ind2] / (pars_true[i,j,3] / null_tau) 
            T_stats[i,j,k] = sqrt(n) * (beta_est[ind] - pars_true[i,j,1] * beta[ind])  / (pars_true[i,j,3] / tau) 
        end 
    end 
end 

# Constructing tables 
quantiles = [0.01 .05 .1 .25 .5 .75 .9 .95 .99]
normal_quantile_vals = quantile.(Normal(), quantiles)
chisq_quantile_vals = quantile.(Chisq(dfs), quantiles)

normal_quantile_table = Matrix{Any}(undef, 2 * length(kappas) * length(intercepts), 3 + length(normal_quantile_vals))
row_count = 1
for (i, kappa) in enumerate(kappas)
    for (j, intercept) in enumerate(intercepts) 
        normal_quantile_table[row_count:row_count+1, 1] = ["null", "non-null"]
        normal_quantile_table[row_count:row_count+1, 2] .= kappa  
        normal_quantile_table[row_count:row_count+1, 3] .= intercept 
        for (k, q) in enumerate(normal_quantile_vals) 
            normal_quantile_table[row_count, k + 3] = mean(null_T_stats[i, j, :] .<= q) 
            normal_quantile_table[row_count + 1, k + 3] = mean(T_stats[i, j, :] .<= q) 
        end 
        row_count +=2 
    end 
end 
nams = ["Coordinate", "κ", "β₀", "%Z <= " .* string.(round.(normal_quantile_vals, digits = 3))...]
normal_quantile_table = DataFrame(normal_quantile_table, nams) 
chisq_quantile_table = Matrix{Any}(undef, length(kappas) * length(intercepts), 2 + length(chisq_quantile_vals))
row_count = 1
for (i, kappa) in enumerate(kappas)
    for (j, intercept) in enumerate(intercepts) 
        chisq_quantile_table[row_count, 1] = kappa  
        chisq_quantile_table[row_count, 2] = intercept 
        for (k, q) in enumerate(chisq_quantile_vals) 
            chisq_quantile_table[row_count, k + 2] = mean(llr[i, j, :] .<= q) 
        end 
        row_count += 1
    end 
end 
nams = ["κ", "β₀", "%LLR <= " .* string.(round.(chisq_quantile_vals, digits = 3))...]
chisq_quantile_table = DataFrame(chisq_quantile_table, nams)  

CSV.write(joinpath(results_path, "intercept-simul-z-stat-quants.csv"), normal_quantile_table)
CSV.write(joinpath(results_path, "intercept_simul-plr-stat-quants.csv"), chisq_quantile_table)

# adjusted Z-statistics table null coordinate
nice_col = [(chisq_quantile_table[i,"β₀"] == 1.5 ? "κ:" * string.(chisq_quantile_table[i,"κ"]) * ", " : "       ") * "β₀:" * string.(chisq_quantile_table[i,"β₀"]) for i in axes(chisq_quantile_table, 1)] 
nice_df = hcat(nice_col, 100 .* filter(row -> row.Coordinate == "null", normal_quantile_table)[:, 4:end])
headers = (vcat("Setting", DataFrames.names(nice_df)[2:end]), vcat("", vec(string.(convert.(Int, quantiles * 100))) .* "% quantile"))
pretty_table(nice_df, body_hlines = [5, 10, 15], header = headers)

# adjusted Z-statistics table non-null
nice_col = [(chisq_quantile_table[i,"β₀"] == 1.5 ? "κ:" * string.(chisq_quantile_table[i,"κ"]) * ", " : "       ") * "β₀:" * string.(chisq_quantile_table[i,"β₀"]) for i in axes(chisq_quantile_table, 1)] 
nice_df = hcat(nice_col, 100 .* filter(row -> row.Coordinate == "non-null", normal_quantile_table)[:, 4:end])
headers = (vcat("Setting", DataFrames.names(nice_df)[2:end]), vcat("", vec(string.(convert.(Int, quantiles * 100))) .* "% quantile"))
pretty_table(nice_df, body_hlines = [5, 10, 15], header = headers)

# adjusted PLR statistics table 
nice_df = hcat(nice_col, 100 .* chisq_quantile_table[:, 3:end])
headers = (vcat("Setting", DataFrames.names(nice_df)[2:end]), vcat("", vec(string.(convert.(Int, quantiles * 100))) .* "% quantile"))
pretty_table(nice_df, body_hlines = [5, 10, 15], header = headers)

