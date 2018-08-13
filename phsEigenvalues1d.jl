using PyPlot
using LinearAlgebra
using Printf

include("phs1.jl")

###########################################################################

function PHSorFDeigs(; ptb=30, frq=3, dx=1/16, dt=dx/8, pol=5,
    showPlots=true)
    
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
    
    alp = (-1+2*mod(K,2)) * 2.0 ^ -(1+3*K)
    # if K == 1
    #     alp =  2.0 ^ -4
    # elseif K == 2
    #     alp = -2.0 ^ -7
    # elseif K == 3
    #     alp =  2.0 ^ -10
    # elseif K == 4
    #     alp = -2.0 ^ -13
    # else
    #     error("Haven't considered this K yet.")
    # end
    
    Wx  = getPeriodicDM(z=x, x=x, m=1,
        phs=phs, pol=pol, stc=stc, period=2)
    
    Whv = getPeriodicDM(z=x, x=x, m=2*K,
        phs=phs, pol=pol, stc=stc, period=2)
    
    W = -Wx + alp*dx^(2*K-1)*Whv
    maxErr = maximum(abs.(W*f(x)+fp(x)))
    
    if showPlots
        
        e1 = eigen(Matrix(-Wx*dt)).values
        maxReal1 = maximum(real(e1))
        
        e2 = eigen(Matrix(W*dt)).values
        maxReal2 = maximum(real(e2))
        
        figure(1)
        subplot(121)
        scatter(real(e1), imag(e1), color="red")
        grid(:true)
        title(@sprintf("maxReal=%g", maxReal1))
        subplot(122)
        scatter(real(e2), imag(e2), color="black")
        grid(:true)
        title(@sprintf("maxReal=%g", maxReal2))
        show()
        
    end
    
    return x, W
    
end

###########################################################################
