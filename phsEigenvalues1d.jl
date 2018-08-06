using Plots
gr()
include("phs1.jl")

###########################################################################

function plotPHSeigs(; frq=3, dx=1/16, pol=5, FD=false)
    
    f(x) = sin.(frq*pi*x)
    fp(x) = frq*pi*cos.(frq*pi*x)
    
    x = -1+dx/2 : dx : 1-dx/2
    s = .30 * dx
    x = x + s*( -1 + 2*rand(size(x)) )
    
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
        alp = -2.0 ^ -8
    elseif K == 3
        alp =  2.0 ^ -12
    elseif K == 4
        alp = -2.0 ^ -16
    else
        error("Haven't considered this K yet.")
    end
    
    Wx  = getPeriodicDM(z=x, x=x, m=1,
        phs=phs, pol=pol, stc=stc, period=2.0)
    Whv = getPeriodicDM(z=x, x=x, m=2*K,
        phs=phs, pol=pol, stc=stc, period=2.0)
    
    W = -Wx + alp*dx^(2*K-1)*Whv
    println(maximum(abs.(W*f(x)+fp(x))))
    
    D = eig(full(W))[1]
    e = diagm(D)
    println(maximum(real(e)))
    
    scatter(real(e), imag(e), legend=false, color="black")

end

###########################################################################
