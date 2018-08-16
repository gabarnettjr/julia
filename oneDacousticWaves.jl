using PyPlot
using Printf
using LinearAlgebra
using SparseArrays

include("rk.jl")
include("phs1.jl")

function main(; plotSolution=false, plotEigenvalues=true,
periodic=false, ptb=20, dx=1/32, phs=7, pol=5, stc=11, alp=2^-5)
    
    tf = 1
    dt = dx / 4
    frq = 2
    rks = 3
    
    x = -1-dx/2 : dx : 1+dx/2
    s = ptb/100 * dx
    x = x + s*(-1 .+ 2*rand(Float64, size(x)))

    xi = x[2:end-1]
    ni = length(xi)

    U0 = [cos.(frq*pi*xi); zeros(size(xi))]

    U = copy(U0)
    
    #######################################################################
    
    #Boundary condition enforcement weights:
    
    wDirichletLeft = getWeights(z=-1, x=x[1:stc],
        m=0, phs=phs, pol=pol)
    wDirichletLeft = -wDirichletLeft[2:end] ./ wDirichletLeft[1]
    
    wDirichletRight = getWeights(z=1,  x=x[end:-1:end-stc+1],
        m=0, phs=phs, pol=pol)
    wDirichletRight = -wDirichletRight[2:end] ./ wDirichletRight[1]
    
    wNeumannLeft = getWeights(z=-1, x=x[1:stc],
        m=1, phs=phs, pol=pol)
    wNeumannLeft = -wNeumannLeft[2:end] ./ wNeumannLeft[1]
    
    wNeumannRight = getWeights(z=1,  x=x[end:-1:end-stc+1],
        m=1, phs=phs, pol=pol)
    wNeumannRight = -wNeumannRight[2:end] ./ wNeumannRight[1]
    
    #######################################################################
    
    #Main weights for approximating spatial derivatives:

    WxPeriodic = Matrix(getPeriodicDM(z=xi, x=xi, m=1,
        phs=phs, pol=pol, stc=stc, period=2))

    WxDirichlet = Matrix(getDM(z=xi, x=x, m=1,
        phs=phs, pol=pol, stc=stc))

    WxNeumann = copy(WxDirichlet)

    WhvPeriodic = Matrix(getPeriodicDM(z=xi, x=xi, m=pol+1,
        phs=phs, pol=pol, stc=stc, period=2))

    WhvDirichlet = Matrix(getDM(z=xi, x=x, m=pol+1,
        phs=phs, pol=pol, stc=stc))

    WhvNeumann = copy(WhvDirichlet)

    for j in 2 : stc
        for i in 1 : ni
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
    WhvPeriodic = alp*dx^(pol+1) * WhvPeriodic
    WhvDirichlet = alp*dx^(pol+1) * WhvDirichlet[:,2:end-1]
    WhvNeumann = alp*dx^(pol+1) * WhvNeumann[:,2:end-1]
    
    #######################################################################
    # #TEST WxDirichlet and WxNeumann:
    # figure()
    # k = 3
    # plot(xi, WxDirichlet*sin.(k*pi*xi) - k*pi*cos.(k*pi*xi))
    # plot(xi, WxNeumann*cos.(k*pi*xi) + k*pi*sin.(k*pi*xi))
    # # plot(xi, fp(xi))
    # show()
    #######################################################################

    #Get block matrix and possibly plot eigenvalues:

    A = zeros(2*ni, 2*ni)
    Aperiodic = zeros(2*ni, 2*ni)

    Aperiodic[1:ni,1:ni] = WhvPeriodic
    Aperiodic[1:ni, ni+1:end] = WxPeriodic
    Aperiodic[ni+1:end, 1:ni] = WxPeriodic
    Aperiodic[ni+1:end, ni+1:end] = WhvPeriodic

    A[1:ni, 1:ni] = WhvNeumann
    A[1:ni, ni+1:end] = WxDirichlet
    A[ni+1:end, 1:ni] = WxNeumann
    A[ni+1:end, ni+1:end] = WhvDirichlet

    if plotEigenvalues
        ei = eigen(dt*A).values
        eiPeriodic = eigen(dt*Aperiodic).values
        figure(2)
        clf()
        subplot(121)
        plot(real(ei), imag(ei), marker=".", color="black",
            linestyle="none")
        title(@sprintf("nonperiodic, maxReal=%g", maximum(real(ei))))
        subplot(122)
        plot(real(eiPeriodic), imag(eiPeriodic), marker=".", color="black",
            linestyle="none")
        title(@sprintf("periodic, maxReal=%g", maximum(real(eiPeriodic))))
        show()
    end

    #######################################################################
    
    #Time-stepping stuff:

    t = 0
    nTimesteps = tf / dt

    dUdt = zeros(2*ni)
    function odefun(t, U, dUdt)
        if periodic
            dUdt = Aperiodic * U
        else
            dUdt = A * U
        end
        return dUdt
    end
    
    q1 = dUdt
    q2 = zeros(size(dUdt))
    if rks == 4
        q3 = zeros(size(dUdt))
        q4 = zeros(size(dUdt))
    end

    for i in 1 : nTimesteps
        if rks == 3
            t, U = rk3!(t, U, odefun, dt, q1, q2)
        elseif rks == 4
            t, U = rk4!(t, U, odefun, dt, q1, q2, q3, q4)
        else
            error("Choose rks=3 or rks=4 for this problem.")
        end
        if maximum(abs.(U)) > 2
            error("unstable in time.")
        end
        if plotSolution
            if mod(i, Int64(1/16/dx)) == 0
                figure(3)
                plot(xi, U[1:ni])
                axis(:image)
                axis([-1,1,-1,1])
                title(@sprintf("numerical solution, t=%1.2f",t))
                show()
                pause(.01)
                if i != nTimesteps
                    clf()
                end
            end
        end
    end

    figure(4)
    plot(xi, U[1:ni] - U0[1:ni])
    title("error at final time")
    show()

end
