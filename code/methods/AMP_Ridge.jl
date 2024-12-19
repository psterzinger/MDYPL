## Philipp Sterzinger 23.10.2023
## Provided as is, although state evolutions have been unit tested, no guarantee for correctness is given 

module AMP_Ridge 

export find_params_ridge_nonlinearsolve, find_params_ridge_nlsolve, logistic_ridge

using LinearAlgebra, Statistics, Cuba, NLsolve, SciMLBase, NonlinearSolve, FiniteDiff, Optim, FastGaussQuadrature

include("mDYPL_helpers.jl")
include("logistic_loglikl.jl")

## Approximate state evolution equations
function eq_bin!(res, param, kappa, gamma, lambda; verbose = false,  int_abstol = 1e-10, int_reltol = 1e-10)
    
    mu = param[1]
    b = param[2] 
    sigma = param[3]
    if verbose 
        println("μ = ", mu, ", b = ", b, ", σ = ", sigma)
    end 

    res .= Cuba.cuhre((x,g) -> eq_bin_integrands!(g,x,gamma,kappa,lambda,mu,sigma,b),2,3; atol = int_abstol, rtol = int_reltol).integral
    if verbose
        println("Max residual: $(maximum(abs.(res)))")
    end  
end 

function eq_bin_integrands!(g, x, gamma, kappa, lambda, mu, sigma, b; cutoff = 50000) 
    jacobian = bivar_coord_transform!(x) 
    Z = x[1]
    G = x[2]
    tdens = jacobian * dnorm(Z) * dnorm(G)

    Z = -gamma * Z
    Z_s = sqrt(kappa) * sigma * G - mu * Z

    ∂_ζ = ∂ζ(Z) 
    prox_val = prox_bζ(Z_s, b)
    b∂ζ = b * ∂ζ(prox_val)
    pre_mult = 2.0 *  ∂_ζ * tdens 

    ## mu 
    g[1] = 2.0 * tdens * ∂²ζ(Z) * prox_val + mu * kappa
    ## b 
    g[2] = pre_mult / (1.0 + b * ∂²ζ(prox_val)) - 1.0 + kappa - lambda * b
    ##sigma 
    g[3] = pre_mult * (Z_s - prox_val)^2 - (sigma * kappa)^2
    g[abs.(g).>cutoff] .= NaN 
end 

"""
    find_params_ridge_nlsolve(kappa,gamma,lambda; kwargs...) 

Solve the stationary logistic ridge state evolution equations for `mu, b, sigma` using `NLsolve` given parameters `kappa`, `gamma`, `lambda` and return `NLsolve` return struct. 
N.b.: `sigma` corresponds to `σ̄ / √κ` in Salehi et al. (2019)

# Arguments 
`kappa::Float64`: Asymptotic ratio of columns/rows of design matrix ∈ (0,1)

`gamma::Float64`: Asymptotic variance of true unobserved linear predictors ≥ 0.0

`lambda::Float64`: Shrinkage paramater lambda of ridge penalty ≥ 0.0

# Keyword Arguments 
`verbose::Bool=false`: print solver information at each Newton step

`x_init::Union{Missing,Vector{Float64}}=missing`: provide custom starting values of mu, b, sigma

`constrained_solve::Bool=false`: use constrained solver functionality of `NLsolve`

`reformulation::Symbol=:smooth`: use one of `NLsolve`'s constrained_solve reformulation options `:smooth`, `:minmax`

`method::Symbol=:newton`: provide one of `NLsolve`'s solvers `:newton`, `:trust_region`, `:anderson`

`iterations::Int64=100`:  maximum number of iterations

`linesearch::Any=LineSearches.BackTracking(order=3)`: linesearch algorithm from `LineSearches`

`ftol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`int_abstol::Float64=1e-10` and `int_reltol::Float64 = 1e-10`: the requested relative (ϵ_rel) and absolute (ϵ_abs​) accuracies of the integrals; see `Cuba` for more info 

# Examples 
```jldoctest
julia> sol = find_params_ridge_nlsolve(0.2,sqrt(5),0.0; x_init = [1.5,3.,4.7]);
julia> sol.zero 
3-element Vector{Float64}:
 1.4993501154076798
 3.0269261421855505
 4.743553128667718
```

See also [`find_params_ridge_nonlinearsolve`](@ref)
"""
function find_params_ridge_nlsolve(kappa, gamma, lambda;
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

    if ismissing(x_init) 
        x_init = vcat(1.,1.,1 + sqrt(kappa) * gamma)
    end 

    f_init =  Vector{Float64}(undef,3) 
    lower = [0.0,0.0,0.0] .+ eps()
    upper = [Inf,Inf,Inf]

    # Setup system of equations & Jacobian 
    f!(r,param) = eq_bin!(r, param, kappa, gamma, lambda; verbose = verbose, int_abstol = int_abstol, int_reltol = int_reltol)
    cache = FiniteDiff.JacobianCache(x_init) 
    J!(J,param) = FiniteDiff.finite_difference_jacobian!(J,f!,param,cache)
    df = OnceDifferentiable(f!, J!, x_init, f_init)

    # solve 
    if verbose 
        println("Solve parameters for: κ = $kappa, γ = $gamma, λ = $lambda \n")
    end 
    if constrained_solve
        sol = mcpsolve(df,lower,upper,x_init, reformulation = reformulation, method = method, iterations = iterations, linesearch = linesearch, ftol = ftol)
    else 
        sol = nlsolve(df, x_init, method = method, iterations = iterations, linesearch = linesearch, ftol = ftol) 
    end 

    if sol.zero[end] < 0 & converged(sol) 
        x_init = sol.zero 
        x_init[end] *= -1 
        #println("Got negative standard deviation, try again...")
        if constrained_solve
            sol = mcpsolve(df,lower,upper,x_init, reformulation = reformulation, method = method, linesearch = linesearch, ftol = ftol)
        else 
            sol = nlsolve(df, x_init, method = method, linesearch = linesearch, ftol = ftol) 
        end 
    end 
    sol
end 

"""
    find_params_ridge_nonlinearsolve(kappa,gamma,lambda; <keyword arguments>) 

Solve the stationary mDYPL state evolution equations for `mu, b, sigma` using `NonlinearSolve` given parameters `kappa`, `gamma`, `lambda` and return `Vector{Float64}`. 
N.b.: `sigma` corresponds to `σ̄ / √κ` in Salehi et al. (2019)

# Arguments 
`kappa::Float64`: Asymptotic ratio of columns/rows of design matrix ∈ (0,1)

`gamma::Float64`: Asymptotic variance of true unobserved linear predictors ≥ 0.0

`lambda::Float64`: Shrinkage paramater lambda of ridge penalty ≥ 0.0

# Keyword Arguments 
`verbose::Bool=false`: print solver information at each Newton step

`x_init::Union{Missing,Vector{Float64}}=missing`: provide custom starting values of mu, b, sigma

`method::SciMLBase.AbstractNonlinearAlgorithm=NewtonRaphson()`: provide one of `NonlinearSolve`'s solvers `NewtonRaphson()`, `TrustRegion()`

`iterations::Int64=100`: maximum number of iterations

`abstol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`reltol::Float64=1e-6`: infinite norm of residuals under which convergence is declared

`specialize::DataType=SciMLBase.FullSpecialize`: control the amount of compilation specialization is performed for the NonlinearProblem; see `SciMLBase`

`int_abstol::Float64=1e-10` and `int_reltol::Float64 = 1e-10`: the requested relative (ϵ_rel) and absolute (ϵ_abs​) accuracies of the integrals; see `Cuba` for more info

# Examples 
```jldoctest
julia> sol = find_params_ridge_nonlinearsolve(0.2,sqrt(5),0.0)
3-element Vector{Float64}:
 1.4993500794510768
 3.026926055407923
 4.7435529873394415
```

See also [`find_params_ridge_nlsolve`](@ref)
"""
function find_params_ridge_nonlinearsolve(kappa, gamma, lambda;
    verbose::Bool = false,
    x_init::Union{Missing,Vector{Float64}} = missing, 
    method::SciMLBase.AbstractNonlinearAlgorithm = NewtonRaphson(), 
    iterations::Int64 = 100,
    abstol::Float64 = 1e-6,
    reltol::Float64 = 1e-6, 
    specialize::DataType = SciMLBase.FullSpecialize,
    int_abstol::Float64 = 1e-10, 
    int_reltol::Float64 = 1e-10, 
    )::Vector{Float64}

    if ismissing(x_init) 
        x_init = vcat(1.,1., 1. + sqrt(kappa) * gamma)  
    end 

    # Setup system of equations & Jacobian 
    jp = Matrix{Float64}(undef,3,3) 

    f!(r,param,lambda) = eq_bin!(r, param, kappa, gamma, lambda; verbose = verbose::Bool, int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
    cache = FiniteDiff.JacobianCache(x_init) 
    g!(r,param) = eq_bin!(r, param, kappa, gamma, lambda; verbose = verbose::Bool, int_abstol = int_abstol::Float64, int_reltol = int_reltol::Float64)
    J!(J,param,lambda) = FiniteDiff.finite_difference_jacobian!(J,g!,param,cache)

    function Jv!(Jv,vec,param,lambda) 
        J!(J,param,lambda)
        mul!(Jv,J,vec) 
        return nothing 
    end 

    function Jvt!(Jv,vec,param,lambda) 
        J!(J, param,lambda) 
        mul!(Jv,J',vec)  
        return nothing 
    end

    # define nonlinear problem
    NF = NonlinearFunction{true, specialize}(f!; 
        jac = J!, 
        jvp = Jv!, 
        vjp = Jvt!, 
        jac_prototype = jp
    )

    # solve 
    if verbose 
        println("\n Solve parameters for: κ = $kappa, γ = $gamma, λ = $lambda \n")
    end 
    probN = NonlinearProblem(NF, x_init, lambda)
    sol = Vector{Float64}(solve(probN, method, reltol = reltol, abstol = abstol, maxiters = iterations))

    if sol[end] < 0
        x_init = sol
        x_init[end] *= -1 
        probN = NonlinearProblem(NF, x_init, lambda)
        sol = Vector{Float64}(solve(probN, method, reltol = reltol, abstol = abstol, maxiters = iterations))
    end 

    sol
end 

## logistic ridge estimator
"""
    logistic_ridge(y, X, lambda; beta_init = missing, kwargs...) 

Compute the ridge estimator for a logistic regression model using data `y`,`X` and shrinkage parameter `alpha` and return a `Optim.optimize` return struct. 

# Arguments 
`y::Vector`: Vector of binary responses

`X::Matrix`: Matrix of covariates 

`lambda::Float64`: Shrinkage paramater lambda of ridge penalty ≥ 0.0

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
julia> lambda = 0.1;  # Regularization parameter
julia> lr = logistic_ridge(y, X, lambda; beta_init = beta);  # Fit the model
julia> Optim.minimizer(lr)
100-element Vector{Float64}:
 0.007007540329276629
 ⋮
 0.2993768841574827
```
"""
function logistic_ridge(y, X, lambda; beta_init = missing, kwargs...) 
    (n,p) = size(X)
    mu_buff = Vector{Float64}(undef,n)
    eta = similar(mu_buff) 
    X_buff = similar(X) 

    f(beta) = ridge_loglikl(beta, y, X, lambda, n, p, eta, mu_buff)  
    g!(g,beta) = ridge_loglikl_grad!(g, beta, y, X, lambda, n, p, eta, mu_buff)  
    h!(H,beta) = ridge_loglikl_hess!(H, beta, X, lambda, n, p, eta, mu_buff, X_buff)  

    if ismissing(beta_init)
        beta_init = zeros(size(X,2)) 
    end 
    Optim.optimize(f, g!, h!, beta_init; kwargs...) 
end 

end 

