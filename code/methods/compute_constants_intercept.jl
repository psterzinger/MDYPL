using RCall, Optim, LinearAlgebra, Statistics, NonlinearSolve, StatsBase;

intercept_unknown = true


include(joinpath("/Users/yiannis/Repositories/hdl/code/mDYPL_supplementary", "Scripts", "mDYPL.jl"))
using .mDYPL

enp_file = ARGS[1]

@rput enp_file
R"""
load(enp_file)
"""
@rget aux start warm

nsettings = size(aux)[1]
consts = Matrix{Float64}(undef, nsettings, 4)
funs = Matrix{Float64}(undef, nsettings, 4)
o = ones(4)

macro threads_if(flag, expr)
    quote
        if $(flag)
            Threads.@threads $expr
        else
            $expr
        end
    end |> esc
end

@threads_if !warm for i in 1:nsettings
    kappa = aux[i, "kappa"]
    eta = aux[i, "eta"]
    alpha = aux[i, "alpha"]
    beta0hat = aux[i, "beta0hat"]

    cstart = start

    k = 0
    for j in 1:size(cstart)[2]
        k = k + 1
        consts[i, :] = mDYPL.solve_mDYPL_SE(kappa, alpha, eta, beta0hat,
                                            use_eta = true,
                                            intercept_unknown = intercept_unknown,
                                            method = TrustRegion(),
                                            abstol = 1e-10,
                                            reltol = 1e-10,
                                            verbose = false,
                                            x_init = cstart[:, j])
        mDYPL.mDYPL_SE!(o, consts[i, :], kappa, eta, alpha, beta0hat, use_eta = true, intercept_unknown = intercept_unknown,
                            with_intercept = true)

        if sum(abs.(o)) < 1e-06
            break
        end
            
    end

    if warm & i > 1
        cstart = start[:, 1]
        cstart[:, 1] = consts[i, :]
    end
    
   
    if consts[i, 3] < 0 
        consts[i, :] = mDYPL.solve_mDYPL_SE(kappa, alpha, eta, beta0hat,
                                            use_eta = true,
                                            intercept_unknown = intercept_unknown,
                                            method = TrustRegion(),
                                            abstol = 1e-10,
                                            reltol = 1e-10,
                                            verbose = false,
                                            x_init = vcat(consts[i, 1:2], -consts[i, 3], consts[i, 4]))
    end

    funs[i, :] = o
    o[:] = ones(4)
    println(i, 
            " | κ = ", round(kappa, digits = 3),
            " | α = ", round(alpha, digits = 3),
            " | η = ", round(eta, digits = 3),
            " | ", round.(consts[i, :], digits = 3),
            " | ", round.(funs[i, :], digits = 10),
            " | init = ", cstart[:, k])
end

@rput consts funs enp_file
R"""
consts <- data.frame(consts)
funs <- data.frame(funs)
names(consts) <- c("mu", "b", "sigma", "intercept")
names(funs) <- c("mu_fun", "b_fun", "sigma_fun", "intercept_fun")
save(consts, funs, file = enp_file)
"""
