n_cores = 10

using Random, Plots, LaTeXStrings, JLD2, ColorSchemes, DataFrames, Distributed, Plots.Measures, RCall

addprocs(n_cores) 

@everywhere begin 
    using LinearAlgebra, Optim, SharedArrays, Random, Statistics, ProgressMeter
    
    supp_path = "." 
    results_path = supp_path * "results" 
    include(joinpath(supp_path, "code", "methods", "mDYPL.jl"))
    @eval using .mDYPL

    function fill_beta_SC(p)
        pp = ceil(Int64, p / 8) 
        vcat(fill(-10,pp),fill(10,pp),fill(0.0,p-2*pp))
    end 

    kappas = (0.1, 0.2, 0.4, 0.6, 0.8, 0.9) 
    gammas = (2, 6, 10, 14, 18)
    #alphas = (0.9, 0.95, 0.99, 0.999) 
    alphas = ["adaptive"] # For a = 1 / (1 + kappa)
    ns = collect(200 : 200 : 8000)
    reps = 1

end 

full_grid = [(0.0, 0.0) for i in  Iterators.product(gammas, kappas)]
for i in eachindex(gammas)
    gamma = gammas[i]
    for j in eachindex(kappas)
        kappa = kappas[j]
        full_grid[i,j] = (kappa, gamma) 
    end 
end 
grid_inds = vcat(1,11,21,7,17,3,13,23,9,19,5,15,25)
grid_vals = [(x,y) = full_grid[ind] for ind in grid_inds]
grid_vals = vcat((0.2,sqrt(0.9)), (sqrt(.9)/22.4, 0.2*22.4), grid_vals[2:end])

betas = [Float64[] for i in eachindex(grid_vals)] 
beta_true = [Float64[] for i in eachindex(grid_vals)] 
beta_length = fill(NaN,length(ns)) 
norms = SharedVector{Float64}(reps) 
params = [Float64[] for i in eachindex(grid_vals)] 
counter = 1 

for a in alphas 
    for i in eachindex(grid_vals)
        (kappa, gamma) = grid_vals[i]
        for k in eachindex(ns) 
            n = ns[k] 
            p = floor(Int64, n * kappa)
            if a == "adaptive"
                a = 1 / (1 + kappa) 
            end 
            beta = fill_beta_SC(p) * gamma / sqrt(kappa) * 0.2 
            beta_true[i] = beta
            println("Setting: $counter / $(length(alphas) * length(grid_vals) * length(ns)), κ: $kappa, γ: $gamma, α: $a, n: $n")
            @showprogress @distributed for l in 1:reps 
                Random.seed!(k * reps + l)
                X = randn(n,p) / sqrt(n) 
                mu = 1.0 ./ (1.0 .+ exp.(.-X*beta)) 
                y = rand(n) .< mu 
                mDYPL = Optim.minimizer(logistic_mDYPL(y,X,a; beta_init = beta)) 
                norms[l] = sum(mDYPL.^2)
            end
            beta_length[k] = mean(norms)  
            counter += 1 
        end 
        betas[i] = copy(beta_length)
    end 
    if length(alphas) > 1
        @save joinpath(results_path, "beta_l2_length_a$a.jld2") betas 
    else 
        @save joinpath(results_path, "beta_l2_length.jld2") betas 
    end 
end 

@rput n_cores supp_path results_path
R"""
library("parallel")
source(file.path(supp_path, "code/methods/compute-pt.R"))

ns <- 200000
set.seed(123)
xzu <- data.frame(X = rnorm(ns), Z = rnorm(ns), U = runif(ns))
pts <- compute_pt(gamma_grid = seq(from = 0.025, to = 22.5, legnth.out = 200), ncores = n_cores, XZU = xzu)
"""
@rget pts 
pts = vcat([1.0, 0.0]', Matrix(pts))

gamma_grid = vcat(0.0,0.001,0.01,collect(0.0:20/100:22.5)[2:end])

scalefontsizes(1.5)
plots = [] 
for a in alphas 
    if a == "adaptive"
        betas = load(joinpath(results_path, "beta_l2_length.jld2"))["betas"] 
    else 
        betas = load(joinpath(results_path, "beta_l2_length_a$a.jld2"))["betas"] 
    end 
    if length(alphas) > 1 
        title = LaTeXString("\$ \\alpha = $(a) \$")
    else 
        title = LaTeXString("\$ \\alpha = 1 / (1 + \\kappa) \$")
    end 
    p = plot(pts[:,1], pts[:,2], 
    label = :none, 
    xlabel = L"$\kappa$", 
    ylabel = L"$\gamma$", 
    linewidth = 0.0, 
    grid=false,
    widen = false, 
    xlims = (0,1), 
    tick_direction = :out, 
    foreground_color_axis = nothing,
    titlefontsize = 11, 
    title = title
    )

    plot!(p, pts[:,1], pts[:,2], 
    fillrange = zero(pts[:,1]),
    fc = :white,  
    lc = :white, 
    label = :none, 
    xlabel = L"$\kappa$", 
    ylabel = L"$\gamma$", 
    )

    plot!(p, pts[:,1], pts[:,2], 
    fillrange = fill(maximum(gamma_grid),length(pts[:,1])),
    fillalpha = 0.15,
    fc = :gray,  
    label = :none, 
    xlabel = L"$\kappa$", 
    ylabel = L"$\gamma$", 
    linewidth = 0.0
    )

    xticks!(p,0:.1:1,vcat(L"0", "",L"0.2","",L"0.4","",L"0.6","",L"0.8","", L"1.0"), z = 1)
    yticks!(p,0:5:20,vcat(L"0",L"5",L"10",L"15",L"20"))

    scatter!(p,grid_vals, color = :gray, label = :none, fillstyle = :none, markersize = 1,)
    xlim = Plots.xlims(p) 
    ylim = Plots.ylims(p) 

    xticklabels = [L"$0$",L"$4\mathrm{k}$",L"$8\mathrm{k}$"]

    for i in eachindex(grid_vals)
        (x,y) = grid_vals[i]
        x_norm = (x  - xlim[1]) / (xlim[2] - xlim[1]) + .005
        y_norm = (y  - ylim[1]) / (ylim[2] - ylim[1]) + .005
        
        y_values = sqrt.(betas[i] ./ collect(200:200:8000))
        
        y_values = sqrt.(betas[i] ./ collect(200:200:8000))
        y_min = minimum(y_values)
        y_max = maximum(y_values)
        y_quarter = y_min + 0.25 * (y_max - y_min)
        y_three_quarter = y_min + 0.75 * (y_max - y_min)

        y_quarter_label = LaTeXString("\$$(round(y_quarter, digits=2))\$")
        y_three_quarter_label = LaTeXString("\$$(round(y_three_quarter, digits=2))\$")

        Plots.plot!(p,
            200:200:8000,
            sqrt.(betas[i] ./ collect(200:200:8000)),
            color = :black,
            fillalpha = .25, 
            alpha = 1, 
            label = :none, 
            markersize = 2, 
            inset = (1, bbox(x_norm, y_norm, 0.15, 0.15, :bottom)),
            subplot = i+1, 
            bg_inside = nothing, 
            #yticks = nothing, 
            xtickfont = font(6), 
            ytickfont = font(6), 
            yticks = ((y_quarter, y_three_quarter), [y_quarter_label, y_three_quarter_label]),
            xticks = ((0,4000,8000), xticklabels)
            #xrotation = 45
            #ticks = nothing
        )
    end
    push!(plots,p)
end 
if length(plots) > 1 
    final_plot = plot(plots...,layout=(2,2), size =  (1200,800), left_margin = 5mm)
    Plots.savefig(final_plot,joinpath(supp_path, "figures", "beta_norm_alphas.pdf"))
else
    final_plot = plots[1] 
    Plots.savefig(final_plot,joinpath(supp_path, "figures", "beta_norm.pdf"))
end 

