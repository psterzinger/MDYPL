## Philipp Sterzinger & Ioannis Kosmidis 07.08.2024
## Provided as is, no guarantee for correctness is given 
module mDYPL

export solve_mDYPL_SE, get_eta, get_gamma, logistic_mDYPL, get_iota

using LinearAlgebra, Statistics, Cuba, NLsolve, SciMLBase, NonlinearSolve, FiniteDiff, Optim, Roots, FastGaussQuadrature

include("mDYPL_helpers.jl")
include("logistic_loglikl.jl")

function eg(gamma, kappa, mu, sigma)
    sqrt(kappa * sigma^2 + mu^2 * gamma^2)
end

function ge(eta, kappa, mu, sigma; min_gamma = 1e-7)
    gamma_sq = (eta^2 - kappa * sigma^2) / mu^2
    sqrt(max(gamma_sq, min_gamma^2))
end

## Approximate state evolution equations
function mDYPL_SE!(res, param, kappa, eta, alpha, intercept = 0; verbose = false,  int_abstol = 1e-10, int_reltol = 1e-10, with_intercept = false, use_eta = true, intercept_unknown = false, success_probability = missing, mu_fixed = false, mu_target = 1.0)
    
    if mu_fixed 
        mu = mu_target 
    else 
        mu = param[1]
    end 

    b = param[2] 
    sigma = param[3]

    if use_eta
        gamma = ge(eta, kappa, mu, sigma)
    else
        gamma = eta 
    end 
    
    if with_intercept
        if intercept_unknown 
            if !use_eta 
                gamma = param[4] 
                b_0 = param[5] 
                if verbose 
                    println("μ = ", mu, ", b = ", b, ", σ = ", sigma, ", γ = ", gamma, ", β₀ = ", b_0) 
                end
                res[1:end-1] .= Cuba.cuhre((x,g) -> mDYPL_SE_integrands!(g,x,gamma,kappa,b_0,alpha,mu,sigma,b,intercept),2,4; atol = int_abstol, rtol = int_reltol).integral
                res[end] = logistic_prob_equation(b_0, gamma, success_probability, ghqr)
            else
                b_0 = param[4] 
                if verbose 
                    println("μ = ", mu, ", b = ", b, ", σ = ", sigma, ", β₀ = ", b_0) 
                end 
                res .= Cuba.cuhre((x,g) -> mDYPL_SE_integrands!(g,x,gamma,kappa,b_0,alpha,mu,sigma,b,intercept),2,4; atol = int_abstol, rtol = int_reltol).integral
            end
        else
            iota = param[4] 
            if verbose
                println("μ = ", mu, ", b = ", b, ", σ = ", sigma, ", ι = ", iota) 
            end 
            res .= Cuba.cuhre((x,g) -> mDYPL_SE_integrands!(g,x,gamma,kappa,intercept,alpha,mu,sigma,b,iota),2,4; atol = int_abstol, rtol = int_reltol).integral
        end 
    else
        if isa(alpha, Function) 
            a = alpha(kappa, gamma) 
        elseif mu_fixed 
            a = param[1]
        else
            a = alpha 
        end 
        if verbose
            if mu_fixed 
                println("α = ", a, ", b = ", b, ", σ = ", sigma)
            else
                println("μ = ", mu, ", b = ", b, ", σ = ", sigma)
            end 
        end 
        res .= Cuba.cuhre((x,g) -> mDYPL_SE_integrands!(g,x,gamma,kappa,a,mu,sigma,b),2,3; atol = int_abstol, rtol = int_reltol).integral
    end
    
    if verbose
        println("Max residual: $(maximum(abs.(res)))")
    end  
end 

function mDYPL_SE_integrands!(g, x, gamma, kappa, a, mu, sigma, b; cutoff = 50000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]
    tdens = jacobian * dnorm(Z) * dnorm(G)

    Z = gamma * Z
    Z_s = mu * Z + sigma * sqrt(kappa) * G

    ∂_ζ = ∂ζ(Z) 
    a_frac = (1 + a) / 2 
    prox_val = prox_bζ(Z_s + a_frac * b, b)
    ∂²_ζ = ∂²ζ(prox_val)
    ∂moreau = a_frac - ∂ζ(prox_val)
    pre_mult = 2.0 * ∂_ζ * tdens 

    ## mu 
    g[1] = pre_mult * Z * ∂moreau
    ## b 
    g[2] = pre_mult / (1.0 + b * ∂²_ζ) - 1.0 + kappa 
    ##sigma 
    g[3] = pre_mult * ∂moreau^2 * b^2 / kappa^2 - sigma^2 

    #g[abs.(g).>cutoff] .= NaN 
end 

function mDYPL_SE_integrands!(g,x,gamma,kappa,intercept,a,mu,sigma,b,iota; cutoff = 50000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]
    tdens = jacobian * dnorm(Z) * dnorm(G)

    Z_s = gamma * mu * Z + sigma * sqrt(kappa) * G + iota 
    Z = gamma * Z + intercept

    ∂_ζ = tdens * ∂ζ(Z) 
    ∂_ζ_n = tdens * ∂ζ(-Z)
    a_frac = (1 + a) / 2 
    prox_val = prox_bζ(Z_s + a_frac * b, b)
    prox_val_n = prox_bζ(a_frac * b - Z_s, b)
    ∂_ζ_prox = ∂ζ(prox_val) 
    ∂_ζ_prox_n = ∂ζ(prox_val_n)  


    ## mu 
    g[1] = ∂_ζ * Z * (a_frac - ∂_ζ_prox) - ∂_ζ_n * Z * (a_frac - ∂_ζ_prox_n)
    ## b 
    g[2] = ∂_ζ / (1 + b * ∂²ζ(prox_val)) + ∂_ζ_n / (1 + b * ∂²ζ(prox_val_n)) - 1 + kappa
    ##sigma 
    g[3] = (∂_ζ * (a_frac - ∂_ζ_prox)^2 + ∂_ζ_n * (a_frac - ∂_ζ_prox_n)^2) * b^2 / kappa^2 - sigma^2  
    ## iota
    g[4] = ∂_ζ * (a_frac - ∂_ζ_prox) - ∂_ζ_n * (a_frac - ∂_ζ_prox_n)
    #g[abs.(g).>cutoff] .= NaN 
end 

function mDYPL_SE_jac_integrands!(vec,x,gamma,kappa,a,mu,sigma,b; cutoff = 50000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]
    tdens = jacobian * dnorm(Z) * dnorm(G)

    Z = gamma * Z
    Z_s = mu * Z + sigma * sqrt(kappa) * G

    ∂_ζ = ∂ζ(Z) 
    a_frac = (1 + a) / 2 
    prox_val = prox_bζ(Z_s + a_frac * b, b)
    ∂²_ζ = ∂²ζ(prox_val)
    ζ_frac = ∂²_ζ / (1.0 + b * ∂²_ζ)
    ∂moreau = a_frac - ∂ζ(prox_val)

    ## mu 
    pre_mult = -2.0 * ∂_ζ * Z * ζ_frac * tdens 

    vec[1] = pre_mult * Z #wrt mu
    vec[2] = pre_mult * ∂moreau #wrt b
    vec[3] = pre_mult * sqrt(kappa) * G #wrt sigma

    ## b
    inter_val = b * ∂³ζ(prox_val) / (1.0 + b * ∂²_ζ)
    pre_mult = -2.0 * ∂_ζ / (1.0 + b * ∂²_ζ)^2 * tdens

    vec[4] = pre_mult * Z * inter_val #wrt mu
    vec[5] = pre_mult * (∂²_ζ + ∂moreau * inter_val) #wrt b
    vec[6] = pre_mult * sqrt(kappa) * G *inter_val #wrt sigma
    
    ### sigma 
    pre_mult = -4.0 * b^2 / kappa^2 * ∂_ζ * ∂moreau * ζ_frac * tdens
    
    vec[7] = pre_mult * Z #wrt mu
    vec[8] = pre_mult * ∂moreau #wrt b
    vec[9] = pre_mult * sqrt(kappa) * G #wrt sigma

    vec[8] += 4 * b / kappa^2 * ∂_ζ * ∂moreau^2 * tdens 
    vec[9] += -2.0 * sigma

   #vec[abs.(vec).>cutoff] .= NaN 
   return nothing 
end 

function mDYPL_SE_jac!(J, param, kappa, gamma, a, vec; int_abstol = 1e-10, int_reltol = 1e-10) 
    mu = param[1]
    b = param[2] 
    sigma = param[3]
   
    vec.= Cuba.cuhre((x,g) -> mDYPL_SE_jac_integrands!(g, x, gamma, kappa, a, mu, sigma, b),2,9; atol = int_abstol, rtol = int_reltol).integral 
    rowfill!(J,vec) 
end 

function mDYPL_SE_jac_integrands_eta!(vec, x, eta, kappa, a, mu, sigma, b; cutoff = 50000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]
    tdens = jacobian * dnorm(Z) * dnorm(G)

    gamma = ge(eta, kappa, mu, sigma) 

    Z_1 = gamma * Z
    Z_2 = mu * Z_1 + sigma * sqrt(kappa) * G

    ∂_ζ = ∂ζ(Z_1) 
    a_frac = (1 + a) / 2 
    prox_val = prox_bζ(Z_2 + a_frac * b, b)
    ∂²_ζ = ∂²ζ(prox_val)
    ζ_frac = ∂²_ζ / (1.0 + b * ∂²_ζ)
    ∂moreau = a_frac - ∂ζ(prox_val)

    d_sigma =  - (kappa * sigma) / (mu * gamma) * Z
    ∂Z_1_mu = - gamma / mu * Z  
    ∂Z_1_sigma = d_sigma / mu 
    ∂Z_2_sigma = d_sigma + G * sqrt(kappa) 

    partial_Z_1_mu_expr = 2.0 * tdens * (∂_ζ + Z_1 * ∂²ζ(Z_1))
    Z_1_pre_mult = 2.0 * tdens * ∂²ζ(Z_1)

    ## mu 
    pre_mult = -2.0 * ∂_ζ * Z_1 * ζ_frac * tdens 

    vec[1] = ∂moreau * partial_Z_1_mu_expr * ∂Z_1_mu #wrt mu
    vec[2] = pre_mult * ∂moreau #wrt b
    vec[3] = ∂moreau * partial_Z_1_mu_expr * ∂Z_1_sigma + pre_mult * ∂Z_2_sigma #wrt sigma

    ## b
    inter_val = b * ∂³ζ(prox_val) / (1.0 + b * ∂²_ζ)
    pre_mult = -2.0 * ∂_ζ / (1.0 + b * ∂²_ζ)^2 * tdens

    vec[4] = Z_1_pre_mult * ∂Z_1_mu / (1.0 + b * ∂²_ζ) #wrt mu
    vec[5] = pre_mult * (∂²_ζ + ∂moreau * inter_val) #wrt b
    vec[6] = Z_1_pre_mult * ∂Z_1_sigma / (1.0 + b * ∂²_ζ) + pre_mult * inter_val * ∂Z_2_sigma #wrt sigma
    
    ### sigma 
    pre_mult = -4.0 * b^2 / kappa^2 * ∂_ζ * ∂moreau * ζ_frac * tdens
    
    vec[7] = Z_1_pre_mult * ∂Z_1_mu * ∂moreau^2 * b^2 / kappa^2 #wrt mu
    vec[8] = pre_mult * ∂moreau #wrt b
    vec[9] = pre_mult * ∂Z_2_sigma + Z_1_pre_mult * ∂Z_1_sigma * ∂moreau^2 * b^2 / kappa^2#wrt sigma

    vec[8] += 4 * b / kappa^2 * ∂_ζ * ∂moreau^2 * tdens 
    vec[9] += -2.0 * sigma

   #vec[abs.(vec).>cutoff] .= NaN 
   return nothing 
end 

function mDYPL_SE_jac_eta!(J, param, kappa, eta, a, vec; int_abstol = 1e-10, int_reltol = 1e-10) 
    mu = param[1]
    b = param[2] 
    sigma = param[3]

    vec .= Cuba.cuhre((x, g) -> mDYPL_SE_jac_integrands_eta!(g, x, eta, kappa, a, mu, sigma, b), 2, 9; atol = int_abstol, rtol = int_reltol).integral 
    rowfill!(J,vec) 
end 

function mDYPL_SE_jac_mu_fixed!(J, param, kappa, gamma, vec, mu_target; int_abstol = 1e-10, int_reltol = 1e-10)
    a = param[1]
    b = param[2] 
    sigma = param[3]
    vec .= Cuba.cuhre((x,g) -> mDYPL_SE_jac_mu_fixed_vec!(g, x, gamma, kappa, a, mu_target, b, sigma),2 ,9; atol = int_abstol, rtol = int_reltol).integral 
    rowfill!(J,vec) 
end 

function mDYPL_SE_jac_mu_fixed_vec!(vec, x, gamma, kappa, a, mu, b, sigma; cutoff = 5000000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]

    tdens = jacobian * dnorm(Z) * dnorm(G)

    Z = gamma * Z
    Z_s = mu * Z + sigma * sqrt(kappa) * G

    ∂_ζ = ∂ζ(Z) 
    a_frac = (1 + a) / 2 
    prox_val = prox_bζ(Z_s + a_frac * b, b)
    ∂²_ζ = ∂²ζ(prox_val)
    ζ_frac = ∂²_ζ / (1.0 + b * ∂²_ζ)
    ∂moreau = a_frac - ∂ζ(prox_val)

    ## mu 
    pre_mult = -2.0 * ∂_ζ * Z * ζ_frac * tdens 

    vec[1] = Z * ∂_ζ / (1.0 + b * ∂²_ζ) * tdens  #wrt a
    vec[2] = pre_mult * ∂moreau #wrt b
    vec[3] = pre_mult * sqrt(kappa) * G #wrt sigma

    ## b
    inter_val = 0.5 * b^2 * ∂³ζ(prox_val) / (1.0 + b * ∂²_ζ)
    pre_mult = -2.0 * ∂_ζ / (1.0 + b * ∂²_ζ)^2 * tdens

    vec[4] = inter_val * pre_mult #wrt a
    vec[5] = pre_mult * (∂²_ζ + ∂moreau * inter_val) #wrt b
    vec[6] = pre_mult * sqrt(kappa) * G *inter_val #wrt sigma
    
    ### sigma 
    pre_mult = -4.0 * b^2 / kappa^2 * ∂_ζ * ∂moreau * ζ_frac * tdens
    
    vec[7] = 2 * b^2 / kappa^2 * ∂_ζ * ∂moreau / (1.0 + b * ∂²_ζ) * tdens #wrt mu
    vec[8] = pre_mult * ∂moreau #wrt b
    vec[9] = pre_mult * sqrt(kappa) * G #wrt sigma

    vec[8] += 4 * b / kappa^2 * ∂_ζ * ∂moreau^2 * tdens 
    vec[9] += -2.0 * sigma

   #vec[abs.(vec).>cutoff] .= NaN 
   return nothing 
end 

"""
    solve_mDYPL_SE(kappa,a,eta_or_gamma=1.0,intercept=1.0; <keyword arguments>) 
 
Solve the mDYPL state evolution equations using `NonlinearSolve` and return `Vector{Float64}`. 

The following maps of input -> output are supported: 
    - (kappa, a, gamma or eta) -> (mu, b, sigma): `intercept = 0.0`, `use_eta = true` for η parameterisation, `use_eta = false` for γ parameterisation and set `eta_or_gamma` to eta or gamma accordingly
    - (kappa, a, gamma or eta, intercept) -> (mu, b, sigma, iota): `use_eta = true` for η parameterisation, `use_eta = false` for γ parameterisation  and set `eta_or_gamma` to eta or gamma accordingly; if intercept = 0.0, iota = 0.0
    - (kappa, a, iota, succes_probability) -> (mu, b, sigma, gamma, intercept): `intercept = 0.0`, `intercept_unknown = true`, `use_eta = false`, `success_probability` set to estimate of  Pr(y = 1| xᵀβ)
    - (kappa, a, eta, iota) -> (mu, b, sigma, intercept): `intercept = 0.0`, `intercept_unknown = true`, `use_eta = true` and set `eta_or_gamma` to eta 

# Arguments 
`kappa::Float64`: Asymptotic ratio of columns/rows of design matrix ∈ (0,1)

`a::Float64`: Shrinkage paramater of mDYPL estimator ∈ (0,1.0]

`eta_or_gamma::Float64=1.0`: Either γ = lim var(Xᵀβ) or η where η² = μ²γ² + κσ²

`intercept::Float64=0.0`: The intercept parameter of β or the mDYPL estimate of the intercept, if intercept is unknown (keyword `intercept_unknown = true`)

# Keyword Arguments 
`use_eta::Bool=true`: Use η or γ parameterisation of state evolution equations 

`intercept_unknown::Bool=false`: The intercept parameter of β is unknown 

`success_probability::Union{Missing, Float64}`: Estimate of success probability Pr(y = 1| xᵀβ), needed if `intercept_unknown = true` and `use_eta = false`

`verbose::Bool=false`: print solver information at each step

`x_init::Union{Missing,Vector{Float64}}=missing`: provide custom starting values of mu, b, sigma

`method::SciMLBase.AbstractNonlinearAlgorithm=TrustRegion()`: Any of `NonlinearSolve`'s solvers

`iterations::Int64=100`: maximum number of iterations

`abstol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`reltol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`specialize::DataType=SciMLBase.FullSpecialize`: control the amount of compilation specialization is performed for the NonlinearProblem; see `SciMLBase`

`int_abstol::Float64=1e-10` and `int_reltol::Float64 = 1e-10`: the requested relative (ϵ_rel) and absolute (ϵ_abs​) accuracies of the integrals; see `Cuba` for more info

`numeric::Bool=true`: Use jacobian based on finite differences (`FiniteDiff`); Analytic jacobian only available in γ parameterisation without intercept 
"""
# function solve_mDYPL_SE(kappa, eta_or_gamma, a, intercept = 0.0;
#     use_eta::Bool = true, 
#     verbose::Bool = false,
#     x_init::Union{Missing,Vector{Float64}} = missing, 
#     method::Union{SciMLBase.AbstractNonlinearAlgorithm} = NewtonRaphson(), 
#     iterations::Int64 = 100,
#     abstol::Float64 = 1e-6,
#     reltol::Float64 = 1e-6, 
#     specialize::DataType = SciMLBase.FullSpecialize,
#     int_abstol::Float64 = 1e-10, 
#     int_reltol::Float64 = 1e-10, 
#     numeric::Bool = true
#     )::Vector{Float64}

#     with_intercept = intercept != 0.0 

#     if verbose 
#         var_print = use_eta ? "η" : "γ" 
#         if with_intercept
#             println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a, β₀ = $intercept \n")
#         else
#             println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a \n")
#         end 
#     end 
    
#     int_dim = with_intercept ? 4 : 3  

#     if ismissing(x_init) 
#         x_init = !with_intercept ? [1.,1.,1 + sqrt(kappa) * eta_or_gamma] : [1.,1.,1 + sqrt(kappa) * eta_or_gamma, intercept / 2]
#     end 

#     # Setup system of equations & Jacobian 
#     jp = Matrix{Float64}(undef,int_dim,int_dim) 

#     f!(r,param,a) = mDYPL_SE!(r, param, kappa, eta_or_gamma, a, intercept; verbose = verbose::Bool, int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64, with_intercept = with_intercept::Bool, use_eta = use_eta::Bool)

#     if !numeric & !with_intercept & !use_eta
#         v = fill(NaN, 12) 

#         J_analytic!(J,param,a) = mDYPL_SE_jac!(J, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
    
#         function Jv_analytic!(Jv,vec,param,a; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
#             mDYPL_SE_jac!(J, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol) 
#             mul!(Jv,J,vec) 
#             return nothing 
#         end 

#         function Jvt_analytic!(Jv,vec,param,a;  int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
#             mDYPL_SE_jac!(J, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol) 
#             mul!(Jv,J',vec)  
#             return nothing 
#         end

#         # define nonlinear problem
#         NF_analytic = NonlinearFunction{true, specialize}(f!; 
#             jac = J_analytic!, 
#             jvp = Jv_analytic!, 
#             vjp = Jvt_analytic!, 
#             jac_prototype = jp
#         )

#         probN_analytic = NonlinearProblem(NF_analytic, x_init, a)
#         sol = Vector{Float64}(solve(probN_analytic, method, reltol = reltol, abstol = abstol, maxiters = iterations))
#     else
#         cache = FiniteDiff.JacobianCache(x_init) 
#         g!(r,param) = mDYPL_SE!(r, param, kappa,  eta_or_gamma, a, intercept; verbose = verbose::Bool, int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64, with_intercept = with_intercept::Bool, use_eta = use_eta::Bool)
#         J!(J,param,a) = FiniteDiff.finite_difference_jacobian!(J,g!,param,cache)
        
    
#         function Jv!(Jv,vec,param,a) 
#             J!(J,param,a)
#             mul!(Jv,J,vec) 
#             return nothing 
#         end 
    
#         function Jvt!(Jv,vec,param,a) 
#             J!(J, param,a) 
#             mul!(Jv,J',vec)  
#             return nothing 
#         end
    
#         # define nonlinear problem
#         NF = NonlinearFunction{true, specialize}(f!; 
#             jac = J!, 
#             jvp = Jv!, 
#             vjp = Jvt!, 
#             jac_prototype = jp
#         )

#         probN = NonlinearProblem(NF, x_init, a)
#         sol = Vector{Float64}(solve(probN, method, reltol = reltol, abstol = abstol, maxiters = iterations))
#     end

#     sol
# end 

function solve_mDYPL_SE(kappa, a, eta_or_gamma = 1.0, intercept = 0.0;
    use_eta::Bool = false, 
    intercept_unknown::Bool = false, 
    mu_fixed::Bool = false, 
    mu_target::Union{Float64, Int64} = 1, 
    success_probability::Union{Float64, Missing} = missing, 
    verbose::Bool = false,
    x_init::Union{Missing,Vector{Float64}} = missing, 
    method::Union{SciMLBase.AbstractNonlinearAlgorithm} = NewtonRaphson(), 
    iterations::Int64 = 100,
    abstol::Float64 = 1e-6,
    reltol::Float64 = 1e-6, 
    specialize::DataType = SciMLBase.FullSpecialize,
    int_abstol::Float64 = 1e-10, 
    int_reltol::Float64 = 1e-10, 
    numeric::Bool = true
    )::Vector{Float64}

    with_intercept = intercept_unknown || intercept != 0.0 

    ## Checks with regard to intercept unknown, e.g. use eta = false

    if verbose 
        if use_eta 
            var_print = "η" 
        else
            var_print = "γ" 
        end   
        if with_intercept
            if intercept_unknown 
                if use_eta 
                    println("Solve parameters for: κ = $kappa, a = $a, ι̂ = $intercept, $var_print = $eta_or_gamma \n")
                else
                    println("Solve parameters for: κ = $kappa, a = $a, ι̂ = $intercept, p̂ = $success_probability \n")
                end
            else
                println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a, β₀ = $intercept \n")
            end 
        else
            if mu_fixed 
                println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, μ = $mu_target \n")
            else
                println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a \n")
            end 
        end 
    end 
    
    if with_intercept
        if intercept_unknown 
            if !use_eta
                int_dim = 5 
                if ismissing(x_init) 
                    x_init = [1.,1.,1 + sqrt(kappa) * 2, 2 * intercept, 3]
                end 
            else
                int_dim = 4 
                if ismissing(x_init) 
                    x_init = [1.,1.,1 + sqrt(kappa) * eta_or_gamma, 2 * intercept]
                end 
            end 
        else
            int_dim = 4
            if ismissing(x_init) 
                x_init = [1.,1.,1 + sqrt(kappa) * eta_or_gamma, intercept / 2]
            end 
        end 
    else 
        int_dim = 3
        if ismissing(x_init) 
            x_init = [1.,1.,1 + sqrt(kappa) * eta_or_gamma]
        end 
    end 


    # Setup system of equations & Jacobian 
    jp = Matrix{Float64}(undef,int_dim,int_dim) 
    if !numeric & !with_intercept
        f_analytic!(r, param, placeholder) = mDYPL_SE!(r, param, kappa, eta_or_gamma, a, 0.0; verbose = verbose::Bool, int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64,  use_eta = use_eta::Bool, with_intercept = false, success_probability = missing, intercept_unknown = false, mu_fixed = mu_fixed, mu_target = mu_target)
        v = fill(NaN, 9) 
        jac_buff = fill(NaN, int_dim, int_dim)
        if !use_eta 
            if mu_fixed 
                J_analytic_mu_fixed!(J, param, placeholder) = mDYPL_SE_jac_mu_fixed!(J, param, kappa, eta_or_gamma, v, mu_target; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
            
                function Jv_analytic_mu_fixed!(Jv, vec, param, placeholder; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
                    mDYPL_SE_jac_mu_fixed!(jac_buff, param, kappa, eta_or_gamma, v, mu_target; int_abstol = int_abstol, int_reltol = int_reltol) 
                    mul!(Jv, jac_buff, vec) 
                    return nothing 
                end 

                function Jvt_analytic_mu_fixed!(Jv, vec, param, placeholder;  int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
                    mDYPL_SE_jac_mu_fixed!(jac_buff, param, kappa, eta_or_gamma, v, mu_target; int_abstol = int_abstol, int_reltol = int_reltol)
                    mul!(Jv, jac_buff', vec)  
                    return nothing 
                end

                # define nonlinear problem
                NF_analytic_mu_fixed = NonlinearFunction{true, specialize}(f_analytic!; 
                    jac = J_analytic_mu_fixed!, 
                    jvp = Jv_analytic_mu_fixed!, 
                    vjp = Jvt_analytic_mu_fixed!, 
                    jac_prototype = jp
                )
                probN_analytic_mu_fixed = NonlinearProblem(NF_analytic_mu_fixed, x_init, 1.0)
                sol = Vector{Float64}(solve(probN_analytic_mu_fixed, method, reltol = reltol, abstol = abstol, maxiters = iterations))
            else 
                J_analytic!(J, param, placeholder) = mDYPL_SE_jac!(J, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
            
                function Jv_analytic!(Jv, vec, param, placeholder; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
                    mDYPL_SE_jac!(jac_buff, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol) 
                    mul!(Jv, jac_buff, vec) 
                    return nothing 
                end 

                function Jvt_analytic!(Jv, vec, param, placeholder;  int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
                    mDYPL_SE_jac!(jac_buff, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol)
                    mul!(Jv, jac_buff', vec)  
                    return nothing 
                end

                # define nonlinear problem
                NF_analytic = NonlinearFunction{true, specialize}(f_analytic!; 
                    jac = J_analytic!, 
                    jvp = Jv_analytic!, 
                    vjp = Jvt_analytic!, 
                    jac_prototype = jp
                )

                probN_analytic = NonlinearProblem(NF_analytic, x_init, 1.0)
                sol = Vector{Float64}(solve(probN_analytic, method, reltol = reltol, abstol = abstol, maxiters = iterations))
            end
        else
            J_analytic_eta!(J, param, placeholder) = mDYPL_SE_jac_eta!(J, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
        
            function Jv_analytic_eta!(Jv, vec, param, placeholder; int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
                mDYPL_SE_jac_eta!(jac_buff, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol) 
                mul!(Jv, jac_buff, vec) 
                return nothing 
            end 

            function Jvt_analytic_eta!(Jv, vec, param, placeholder;  int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64) 
                mDYPL_SE_jac_eta!(jac_buff, param, kappa, eta_or_gamma, a, v; int_abstol = int_abstol, int_reltol = int_reltol) 
                mul!(Jv, jac_buff', vec)  
                return nothing 
            end

            # define nonlinear problem
            NF_analytic_eta = NonlinearFunction{true, specialize}(f_analytic!; 
                jac = J_analytic_eta!, 
                jvp = Jv_analytic_eta!, 
                vjp = Jvt_analytic_eta!, 
                jac_prototype = jp
            )
            probN_analytic_eta = NonlinearProblem(NF_analytic_eta, x_init, 1.0)
            sol = Vector{Float64}(solve(probN_analytic_eta, method, reltol = reltol, abstol = abstol, maxiters = iterations))
        end
    else
        cache = FiniteDiff.JacobianCache(x_init) 
        f!(r, param, a) = mDYPL_SE!(r, param, kappa, eta_or_gamma, a, intercept; 
            verbose = verbose::Bool, 
            int_abstol = int_abstol::Float64, 
            int_reltol = int_reltol::Float64, 
            with_intercept = with_intercept::Bool, 
            use_eta = use_eta::Bool, 
            success_probability = success_probability::Union{Missing, Float64}, 
            intercept_unknown = intercept_unknown::Bool,
            mu_fixed = mu_fixed::Bool, 
            mu_target = mu_target::Union{Float64, Int64} 
        )

        function g!(r, param) 
            mDYPL_SE!(r, param, kappa,  eta_or_gamma, a, intercept; 
                verbose = verbose::Bool, 
                int_abstol = int_abstol::Float64, 
                int_reltol = int_reltol::Float64, 
                with_intercept = with_intercept::Bool, 
                use_eta = use_eta::Bool, 
                success_probability = success_probability::Union{Missing, Float64}, 
                intercept_unknown = intercept_unknown::Bool, 
                mu_fixed = mu_fixed::Bool, 
                mu_target = mu_target::Union{Float64, Int64} 
            )
        end 
        
        function J!(J, param, a) 
            FiniteDiff.finite_difference_jacobian!(J, g!, param, cache)
        end 
        
        jac_buff = fill(NaN, int_dim, int_dim)
        
        function Jv!(Jv, vec, param, a) 
            J!(jac_buff, param, a)
            mul!(Jv, jac_buff, vec) 
            return nothing 
        end 
    
        function Jvt!(Jv, vec, param, a) 
            J!(jac_buff, param,a) 
            mul!(Jv, jac_buff',vec)  
            return nothing 
        end
    
        # define nonlinear problem
        NF = NonlinearFunction{true, specialize}(f!; 
            jac = J!, 
            jvp = Jv!, 
            vjp = Jvt!, 
            jac_prototype = jp
        )

        probN = NonlinearProblem(NF, x_init, a)
        sol = Vector{Float64}(solve(probN, method, reltol = reltol, abstol = abstol, maxiters = iterations))
    end
    sol
end 

## SLOE 
function get_eta(y, X, alpha; beta_init = missing,  kwargs...)
    ys = alpha * y .+ (1 - alpha) / 2
    beta = Optim.minimizer(logistic_mDYPL(y, X, alpha; beta_init = beta_init, kwargs...))
    lp = X * beta
    probs = @. 1 / (1 + exp(- lp))
    probs[probs .== 1] .= 1 - eps() 
    probs[probs .== 0] .= eps() 
    v = probs .* (1 .- probs)
    XW = X .* sqrt.(v)
    Q = Matrix(qr(XW).Q)
    h = sum(Q .* Q, dims = 2)
    S = @. lp - ((ys - probs) * h) / (v * (1 - h))
    sqrt(var(S))
end

function get_gamma(y, X, alpha;
    eta_hat = missing,
    beta_init = missing,
    verbose = false, 
    x_init = missing, 
    method = TrustRegion(), 
    specialize = SciMLBase.FullSpecialize,
    iterations = 100,
    abstol = 1e-6, 
    reltol = 1e-6, 
    int_abstol = 1e-10, 
    int_reltol = 1e-10,
    kwargs...)
    if ismissing(eta_hat) 
        eta_hat = get_eta(y, X, alpha,
                beta_init = beta_init,
                kwargs...)
    end
    dX = size(X)
    kappa_hat = dX[2] / dX[1]
    params = solve_mDYPL_SE(kappa_hat, eta_hat, alpha;
                            verbose = verbose, 
                            x_init = x_init, 
                            method = method, 
                            iterations = iterations,
                            abstol = abstol, 
                            reltol = reltol, 
                            int_abstol = int_abstol, 
                            int_reltol = int_reltol,
                            specialize = specialize, 
                            use_eta = true
                            )
    o = ones(3) .* 10000
    mDYPL_SE!(o, params, kappa_hat, eta_hat, alpha)
    if (any(abs.(o) .> 0.00001)) 
        out = NaN
    else
        out = ge(eta_hat, kappa_hat, params[1], params[3])
    end
    out
end

function get_iota(y, X, alpha; beta_init = missing,  kwargs...)
    ys = alpha * y .+ (1 - alpha) / 2
    beta = Optim.minimizer(logistic_mDYPL(y, X, alpha; beta_init = beta_init, kwargs...))
    lp = X * beta
    probs = @. 1 / (1 + exp(- lp))
    v = probs .* (1 .- probs)
    XW = X .* sqrt.(v)
    Q = Matrix(qr(XW).Q)
    h = sum(Q .* Q, dims = 2)
    S = @. lp - (ys - probs) / v * h / (1 - h)
    mean(S)
end

## mDYPL estimator
"""
    logistic_mDYPL(y, X, alpha; beta_init = missing, kwargs...) 

Compute the mDYPL estimator for a logistic regression model using data `y`,`X` and shrinkage parameter `alpha` and return a `Optim.optimize` return struct. 

# Arguments 
`y::Vector`: Vector of binary responses

`X::Matrix`: Matrix of covariates 

`alpha::Float64`: Shrinkage paramater of mDYPL estimator ∈ (0,1.0]

# Keyword Arguments 
`beta_init::Union{Missing,Vector{Float64}}=missing`: provide starting values for minimization, if missing `beta_init = zeros(size(X, 2))`

`kwargs...`: keyword arguments to be passed to `Optim.optimize` 

# Examples 
```jldoctest
julia> using Random, Optim  # Load necessary packages
julia> Random.seed!(123);  # Seed the random number generator for reproducibility
julia> n = 1000;  # Number of observations
julia> p = 100;   # Number of features
julia> X = randn(n,p) / sqrt(n);  # Generate a random feature matrix
julia> beta = vcat(fill(0.0, ceil(Int64, p / 2)), fill(10.0, p-ceil(Int64, p / 2)));  # True coefficient vector
julia> y = rand(n) .< 1.0 ./ (1.0 .+ exp.(.-X * beta));  # Generate binary response variable
julia> alpha = 0.1;  # Shrinkage parameter
julia> mDYPL = logistic_mDYPL(y, X, alpha; beta_init = beta);  # Fit the model
julia> Optim.minimizer(mDYPL)
100-element Vector{Float64}:
 -0.056723264412529874
  ⋮
  0.2590704292567825
```
"""
function logistic_mDYPL(y, X, alpha; beta_init = missing, kwargs...) 
    y_star = alpha * y .+ (1-alpha) / 2
    n = length(y)
    mu_buff = Vector{Float64}(undef,n)
    eta = similar(mu_buff) 
    X_buff = similar(X) 

    f(beta) = loglikl(beta, y_star, X, eta, mu_buff)  
    g!(g,beta) = loglikl_grad!(g,beta, y_star, X, eta, mu_buff)  
    h!(H,beta) = loglikl_hess!(H,beta, X, eta, mu_buff, X_buff)  

    if ismissing(beta_init)
        beta_init = zeros(size(X,2)) 
    end 
    Optim.optimize(f, g!, h!, beta_init; kwargs...) 
end 

end 

#= 
"""
    solve_mDYPL_SE_nlsolve(kappa,eta_or_gamma,a; kwargs...) 

Solve the SLOE state evolution equations for `mu, b, sigma, iota` using `NLsolve` given parameters `kappa`, `eta` or `gamma`, `a` and return `NLsolve` return struct. 


# Arguments 

`kappa::Float64`: Asymptotic ratio of columns/rows of design matrix ∈ (0,1)

`eta_or_gamma::Float64`: Either γ = lim var(Xᵀβ₀) or η where η² = μ²γ² + κσ²

`a::Float64`: Shrinkage paramater of mDYPL estimator ∈ (0,1.0]

`intercept::Float64=0.0`: Optional intercept parameter for state evolution equations 

# Keyword Arguments 

`use_eta::Bool=true`: Use η or γ parameterisation of state evolution equations 

`verbose::Bool=false`: print solver information at each Newton step

`x_init::Union{Missing,Vector{Float64}}=missing`: provide custom starting values of mu, b, sigma

`constrained_solve::Bool=false`: use constrained solver functionality of `NLsolve`

`reformulation::Symbol=:smooth`: use one of `NLsolve`'s constrained_solve reformulation options `:smooth`, `:minmax`

`method::Symbol=:newton`: provide one of `NLsolve`'s solvers `:newton`, `:trust_region`, `:anderson`

`iterations::Int64=100`:  maximum number of iterations

`linesearch::Any=LineSearches.BackTracking(order=3)`: linesearch algorithm from `LineSearches`

`ftol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`int_abstol::Float64=1e-10` and `int_reltol::Float64 = 1e-10`: the requested relative (ϵ_rel) and absolute (ϵ_abs​) accuracies of the integrals; see `Cuba` for more info 

"""
function solve_mDYPL_SE_nlsolve(kappa, eta_or_gamma, a, intercept = 0.0;
    use_eta = true, 
    verbose = false, 
    x_init = missing, 
    constrained_solve = true, 
    method = :newton, 
    reformulation = :smooth,
    iterations = 100,
    linesearch = LineSearches.BackTracking(order=3),
    ftol = 1e-6, 
    int_abstol = 1e-10, 
    int_reltol = 1e-10
    )

    with_intercept = intercept != 0.0 
    
    int_dim = with_intercept ? 4 : 3  

    if ismissing(x_init) 
        x_init = !with_intercept ? [1.,1.,1 + sqrt(kappa) * eta_or_gamma] : [1.,1.,1 + sqrt(kappa) * eta_or_gamma, intercept / 2]
    end 


    f_init =  Vector{Float64}(undef,int_dim) 
    lower = fill(0.0, int_dim) .+ eps()
    upper = fill(Inf, int_dim) 


    # Setup system of equations & Jacobian 
    f!(r,param) = mDYPL_SE!(r, param, kappa, eta_or_gamma, a, intercept; verbose = verbose, int_abstol = int_abstol, int_reltol = int_reltol, with_intercept = with_intercept, use_eta = use_eta)
    cache = FiniteDiff.JacobianCache(x_init) 
    J!(J,param) = FiniteDiff.finite_difference_jacobian!(J,f!,param,cache)
    df = OnceDifferentiable(f!, J!, x_init, f_init)


    # solve 
    if verbose 
        var_print = use_eta ? "η" : "γ" 
        if with_intercept
            println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a, β₀ = $intercept \n")
        else
            println("Solve parameters for: κ = $kappa, $var_print = $eta_or_gamma, a = $a \n")
        end 
    end 
    if constrained_solve
        sol = mcpsolve(df,lower,upper,x_init, reformulation = reformulation, method = method, iterations = iterations, linesearch = linesearch, ftol = ftol)
    else 
        sol = nlsolve(df, x_init, method = method, iterations = iterations, linesearch = linesearch, ftol = ftol) 
    end 

    if sol.zero[3] < 0 & converged(sol) 
        x_init = sol.zero 
        x_init[3] *= -1 
        #println("Got negative standard deviation, try again...")
        if constrained_solve
            sol = mcpsolve(df,lower,upper,x_init, reformulation = reformulation, method = method, linesearch = linesearch, ftol = ftol)
        else 
            sol = nlsolve(df, x_init, method = method, linesearch = linesearch, ftol = ftol) 
        end 
    end 
    sol
end 
=#