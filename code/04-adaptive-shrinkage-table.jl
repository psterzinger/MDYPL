using Random, Statistics, PrettyTables, Optim

supp_path = "." 
include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
using .mDYPL

function fill_beta_RA(p)
    pp = ceil(Int64, p / 5) 
    vcat(fill(-3,pp),fill(-3/2,pp),fill(0.0,p-4*pp),fill(3/2,pp),fill(3,pp))
end 

function fill_beta_RA_mod(p)
    pp = ceil(Int64, p / 5) 
    vcat(fill(-3,pp),fill(-3/2,pp),fill(0.0,p-4*pp),fill(3/2,pp),fill(8,pp))
end

n = 1000 
reps = 5000 
count = 0
dfs = []

for i in 1:3 
    if i == 1
        kappa = .2
        p = floor(Int64, kappa * n)
        beta = fill_beta_RA(p) 
    elseif i == 2 
        kappa = 0.05 
        gamma = 10 
        p = floor(Int64, kappa * n)
        beta = fill_beta_RA(p) / sqrt(9/2 * kappa) * gamma
    else 
        kappa = 0.2 
        p = floor(Int64, kappa * n)
        beta = fill_beta_RA_mod(p)
    end 
    a = 1 / (1 + kappa) 
    out = Matrix{Float64}(undef,reps,p) 
    Threads.@threads for j in 1:reps 
        Random.seed!(i * count + j)
        X = randn(n,p) / sqrt(n) 
        mu = 1.0 ./ (1.0 .+ exp.(.-X*beta)) 
        y = rand(n) .< mu 
        out[j,:] = Optim.minimizer(logistic_mDYPL(y, X, a; beta_init = beta))
    end 
    out .-= beta'
    pp = convert(Int64, p / 5)
    coord_means = [mean(out[:,(i-1)*pp+1:i*pp]) for i in 1:5]
    coord_rmse = [sqrt(mean(out[:,(i-1)*pp+1:i*pp].^2)) for i in 1:5]
    aggregate_mean = mean(out) 
    aggregate_rmse = sqrt(mean(out.^2)) 
    res = fill(NaN,2,6) 
    res[1,:] = round.(vcat(coord_means,aggregate_mean), digits = 4) 
    res[2,:] = round.(vcat(coord_rmse,aggregate_rmse), digits = 4) 
    unique_betas = round.(unique(beta), digits = 1) 
    colnames = ["β_{0,j} = $(unique_betas[1])", "β_{0,j} = $(unique_betas[2])", "β_{0,j} = $(unique_betas[3])", "β_{0,j} = $(unique_betas[4])", "β_{0,j} = $(unique_betas[5])", "Aggregate"]
    colnames = reshape(colnames, (1,6))
    if i > 1
        dfs = vcat(dfs, vcat(colnames, res))
    else 
        dfs = vcat(colnames, res)
    end 
    count += 1
end 

dfs = hcat(repeat(["κ:0.2,γ:√0.9", "κ:0.05,γ:10", "κ:0.2,γ:√3.1"], inner = 3), repeat(["","Bias","RMSE"], outer = 3), dfs)
pretty_table(dfs, body_hlines = [3,6], header = fill("",8), formatters = ft_round(4))


