using PyPlot
using LinearAlgebra
using Printf

include("phs1.jl")

###########################################################################

function PHSorFDeigs(; ptb=30, frq=3, dx=1/16, pol=5, showPlots=true)
    
    if mod(pol,2) == 0
        FD = true
    else
        FD = false
    end
    
    f(x) = sin.(frq*pi*x)
    fp(x) = frq*pi*cos.(frq*pi*x)
    
    x = -1+dx/2 : dx : 1-dx/2
    # x = range(-1+dx/2, stop=1-dx/2, step=dx)
    s = ptb/100 * dx
    x = x + s*(-1 .+ 2*rand(Float64, size(x)))
    
    if FD & (mod(pol,2) == 0)
        phs = pol+1
        stc = pol+1
        K = Int64(round(pol/2))
    elseif ~FD & (mod(pol,2) == 1)
        phs = 2*pol + 1
        stc = 2*pol + 1
        K = Int64(round((pol+1)/2))
    else
        error("Invalid parameters.  Choose FD=true and pol=even
            or FD=false and pol=odd.")
    end
    
    if K == 1
        alp =  2.0 ^ -4
    elseif K == 2
        alp = -2.0 ^ -7
    elseif K == 3
        alp =  2.0 ^ -10
    elseif K == 4
        alp = -2.0 ^ -13
    else
        error("Haven't considered this K yet.")
    end
    
    Wx  = getPeriodicDM(z=x, x=x, m=1,
        phs=phs, pol=pol, stc=stc, period=2.0)
    
    Whv = getPeriodicDM(z=x, x=x, m=2*K,
        phs=phs, pol=pol, stc=stc, period=2.0)
    
    W = -Wx + alp*dx^(2*K-1)*Whv
    maxErr = maximum(abs.(W*f(x)+fp(x)))
    
    D = eigen(Matrix(-Wx)).values
    e1 = Matrix(Diagonal(D))
    
    D = eigen(Matrix(W)).values
    e2 = Matrix(Diagonal(D))
    maxReal = maximum(real(e2))
    
    if showPlots
        figure()
        subplot(121)
        scatter(real(e1), imag(e1), color="red")
        title(@sprintf("maxErr=%g", maxErr))
        subplot(122)
        scatter(real(e2), imag(e2), color="black")
        title(@sprintf("maxReal=%g", maxReal))
        show()
    else
        return maxErr, maxReal
    end
    
end

###########################################################################
