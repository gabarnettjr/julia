using PyPlot
using Printf
using LinearAlgebra
using SparseArrays

include("rk.jl")
include("phs1.jl")

function main(; dx=1/32, phs=7, pol=3, stc=7, rks=3)
    
    dt = dx / 8
    frq = 2
    
    x = -1-dx/2 : dx : 1+dx/2
    xi = x[2:end-1]
    
    f(x) = cos.(frq*pi*x)
    # fp(x) = -frq*pi*sin.(frq*pi*x)
    
    U0 = [f(xi); zeros(size(xi))]
    
    #######################################################################
    
    #Boundary condition enforcement weights:
    
    wDirichletLeft = getWeights(z=-1, x=x[1:stc],
        m=0, phs=phs, pol=pol)
    wDirichletLeft = -wDirichletLeft[2:end] ./ wDirichletLeft[1]
    
    wDirichletRight = getWeights(z=1,  x=x[end:-1:end-stc+1],
        m=0, phs=phs, pol=pol)
    wDirichletRight = -wDirichletRight[2:end] ./ wDirichletRight[1]
    
    wNeumannLeft = getWeights(z=-1, x=x[1:stc],
        m=1, phs=phs-2, pol=pol)
    wNeumannLeft = -wNeumannLeft[2:end] ./ wNeumannLeft[1]
    
    wNeumannRight = getWeights(z=1,  x=x[end:-1:end-stc+1],
        m=1, phs=phs-2, pol=pol)
    wNeumannRight = -wNeumannRight[2:end] ./ wNeumannRight[1]
    
    #######################################################################
    
    #Main weights for approximating spatial derivatives:
    
    WxDirichlet = Matrix(getDM(z=xi, x=x, m=1, phs=phs, pol=pol, stc=stc))
    WxNeumann = copy(WxDirichlet)
    WhvDirichlet = getDM(z=xi, x=x, m=pol+1, phs=phs, pol=pol, stc=stc)
    WhvNeumann = copy(WhvDirichlet)
    for j in 2 : stc
        for i in 1 : length(xi)
            WxDirichlet[i,j] = WxDirichlet[i,j] +
                wDirichletLeft[j-1] * WxDirichlet[i,1]
            WxDirichlet[i,end-(j-1)] = WxDirichlet[i,end-(j-1)] +
                wDirichletRight[j-1] * WxDirichlet[i,end]
            WxNeumann[i,j] = WxNeumann[i,j] +
                wNeumannLeft[j-1] * WxNeumann[i,1]
            WxNeumann[i,end-(j-1)] = WxNeumann[i,end-(j-1)] +
                wNeumannRight[j-1] * WxNeumann[i,end]
            WhvDirichlet[i,j] = WhvDirichlet[i,j] +
                wDirichletLeft[j-1] * WhvDirichlet[i,1]
            WhvDirichlet[i,end-(j-1)] = WhvDirichlet[i,end-(j-1)] +
                wDirichletRight[j-1] * WhvDirichlet[i,end]
            WhvNeumann[i,j] = WhvNeumann[i,j] +
                wNeumannLeft[j-1] * WhvNeumann[i,1]
            WhvNeumann[i,end-(j-1)] = WhvNeumann[i,end-(j-1)] +
                wNeumannRight[j-1] * WhvNeumann[i,end]
        end
    end
    WxDirichlet = WxDirichlet[:,2:end-1]
    WxNeumann = WxNeumann[:,2:end-1]
    WhvDirichlet = WhvDirichlet[:,2:end-1]
    WhvNeumann = WhvNeumann[:,2:end-1]
    
    #######################################################################
    #TEST WxDirichlet and WxNeumann:
    figure()
    k = 3
    plot(xi, WxDirichlet*sin.(k*pi*xi) - k*pi*cos.(k*pi*xi))
    plot(xi, WxNeumann*cos.(k*pi*xi) + k*pi*sin.(k*pi*xi))
    # plot(xi, fp(xi))
    show()
    #######################################################################
    
end
