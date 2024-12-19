# logistic loglikl 
function loglikl(beta, y, X, eta, mu_buff)  
    mul!(eta,X,beta) 
    mu_buff .= y .* eta 
    eta .= 1.0 .+ exp.(eta)
    eta .= log.(eta) 
    mu_buff .= eta .- mu_buff

    sum(mu_buff)
end 

function loglikl_grad!(r, beta, y, X, eta, mu_buff)  
    mul!(eta,X,beta) 
    mu_buff .= 1.0 ./ (1.0 .+ exp.(.-eta)) 
    mu_buff .= mu_buff .- y 

    mul!(r,X',mu_buff)
end 

function loglikl_hess!(H, beta, X, eta, mu_buff, X_buff)  
    mul!(eta,X,beta) 
    mu_buff .= 1.0 ./ (1.0 .+ exp.(.-eta)) 
    mu_buff .= mu_buff .* (1.0 .- mu_buff) 
    X_buff .= mu_buff .* X 
    mul!(H,X',X_buff) 
end 

# logistic ridge loglikl 
function ridge_loglikl(beta, y, X, lambda, n, p, eta, mu_buff)  
    mul!(eta,X,beta) 
    mu_buff .= y .* eta 
    eta .= 1.0 .+ exp.(eta)
    eta .= log.(eta) 
    mu_buff .= eta .- mu_buff

    scale = (lambda * n) / (2 * p)

    sum(mu_buff) + scale * sum(beta.^2)
end 

function ridge_loglikl_grad!(r, beta, y, X, lambda, n, p, eta, mu_buff)  
    mul!(eta,X,beta) 
    mu_buff .= 1.0 ./ (1.0 .+ exp.(.-eta)) 
    mu_buff .= mu_buff .- y 

    scale = (lambda * n) / p

    mul!(r,X',mu_buff) 

    r .= r .+ scale * beta 
end 

function ridge_loglikl_hess!(H, beta, X, lambda, n, p, eta, mu_buff, X_buff)  
    mul!(eta,X,beta) 
    mu_buff .= 1.0 ./ (1.0 .+ exp.(.-eta)) 
    mu_buff .= mu_buff .* (1.0 .- mu_buff) 
    X_buff .= mu_buff .* X 
    mul!(H,X',X_buff) 

    scale = (lambda * n) / p

    diag_add!(H,scale)
end 

@inline function diag_add!(H, scalar)
    @inbounds for i in axes(H, 1)
        H[i, i] += scalar
    end
end
