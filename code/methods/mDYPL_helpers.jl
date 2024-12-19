function prox_bζ(x,b; maxiter = 100, br_steps = 5, ϵ = 2*eps())::Float64
    if b < ϵ
        return x 
    elseif -x > log(b - ϵ) - log(ϵ) 
        return x 
    elseif x > log(b - ϵ) - log(ϵ) + b 
        return x - b 
    else 
        z = get_x_init(x,b)
        r = prox_foc(z, b, x)
        if abs(r) < ϵ
            return z
        end
        for _ in 1:maxiter
            delta_z = r / ∂prox_foc(z, b) 
            z -= delta_z
            r = prox_foc(z, b, x)
            if abs(r) < ϵ
                return z
            end
        end
        (l,u) = (-10,10)
        f(z) = prox_foc(z,b,x) 
        fu = f(u)
        fl = f(l)
        c = 1
        while fu * fl > 0 && c < br_steps
            u += 10 
            l -= 10 
            fu = prox_foc(u,b,x) 
            fl = f(l) 
            c += 1
        end 
        if c >= br_steps
            return NaN 
        else
            bisection(f,l,u)[1]
        end
    end  
end

function get_x_init(x,b) 
    if x > 0 
        x - b 
    else 
        x 
    end 
end 

function prox_foc(z, b, x)
    z + b / (1.0 + exp(-z)) - x
end 

function ∂prox_foc(z, b) 
    1.0 + b * exp(z) / (1.0 + exp(z))^2
end 

function bisection(f, a_, b_, atol = 2eps(promote_type(typeof(b_),Float64)(b_)); increasing = sign(f(b_))) # found at https://discourse.julialang.org/t//12658
    a_, b_ = minmax(a_, b_)
    c = middle(a_,b_)
    z = f(c) * increasing
    if z > 0 #
        b = c
        a = typeof(b)(a_)
    else
        a = c
        b = typeof(a)(b_)
    end
    while abs(a - b) > atol
        c = middle(a,b)
        if f(c) * increasing > 0 
            b = c
        else
            a = c
        end
    end
    a, b
end

### logistic link 
function ζ(z) log(1.0 + exp(z)) end 

function ∂ζ(z) 1.0 / (1.0 + exp(-z)) end 

function ∂²ζ(z) 
    mu = 1.0 / (1.0 + exp(-z)) 
    mu * (1.0 - mu)
end 

function ∂³ζ(z) 
    mu = 1.0 / (1.0 + exp(-z)) 
    mu * (1.0 - mu) * (1.0 - 2.0 * mu)
end 

### misc 
function dnorm(x)
    exp(-x^2/2) / sqrt(2 * pi)
end 

function rowfill!(J,vec) 
    count = 1
    @inbounds @simd for i in axes(J,1)
        for j in axes(J,2)
            J[i,j] = vec[count] 
            count += 1 
        end 
    end 
end 

function fullfill!(r,J,vec) 
    p = size(J,2)
    r .= vec[1:p] 
    rowfill!(J,vec[(p+1):end]) 
end 

function bivar_coord_transform!(x)
    u = x[1] 
    u = u * (1.0 - u) 

    v = x[2] 
    v = v * (1.0 - v) 

    x[1] = (2.0 * x[1] - 1.0) / u 
    x[2] = (2.0 * x[2] - 1.0) / v 
    
    (1.0 - 2.0 * u) * (1.0 - 2.0 * v) / (u * v)^2
end 

const ghqr = gausshermite(100, normalize = true)

function logistic_prob_equation(b_0, gamma, prob, quad_rule)
    lp = b_0 .+ gamma .* quad_rule[1]
    pr = 1 ./ (1 .+ exp.(-lp))
    prob - sum(quad_rule[2] .* pr)
end 
