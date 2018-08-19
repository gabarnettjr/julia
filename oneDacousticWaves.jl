using PyPlot
using Printf
using LinearAlgebra
using SparseArrays

include("rk.jl")
include("phs1.jl")

function main(; plotSolution=false, plotEigenvalues=true, plotError=false,
frq=2, ptb=30, periodic=false, dx=1/32, phs=7, pol=5, stc=11, alp=2^-9)
    """
    frq is the frequency of the initial condition
    ptb is the percent perturbation from regularly-spaced nodes
    dx is the average spacing of the nodes
    phs is the exponent of the polyharmonic spline radial function
    pol is the highest degree polynomial in the basis
    stc is the stencil size
    alp is the hyperviscosity parameter
    """
    tf = 10
    dt = dx / 2
    rks = 3
    
    x = -1-dx/2 : dx : 1+dx/2
    s = ptb/100 * dx
    x = x + s*(-1 .+ 2*rand(Float64, size(x)))

    xi = x[2:end-1]
    ni = length(xi)

    U0 = [cos.(pi*frq*xi); zeros(size(xi))]

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

    WxPeriodic = getPeriodicDM(z=xi, x=xi, m=1,
        phs=phs, pol=pol, stc=stc, period=2)

    WxDirichlet = getDM(z=xi, x=x, m=1,
        phs=phs, pol=pol, stc=stc)

    WxNeumann = copy(WxDirichlet)

    WhvPeriodic = getPeriodicDM(z=xi, x=xi, m=phs-1,
        phs=phs, pol=pol, stc=stc, period=2)

    WhvDirichlet = getDM(z=xi, x=x, m=phs-1,
        phs=phs, pol=pol, stc=stc)

    WhvNeumann = copy(WhvDirichlet)

    for j in 2 : stc
        WxDirichlet[:,j] = WxDirichlet[:,j] +
            wDirichletLeft[j-1] * WxDirichlet[:,1]
        WxDirichlet[:,end-(j-1)] = WxDirichlet[:,end-(j-1)] +
            wDirichletRight[j-1] * WxDirichlet[:,end]
        WxNeumann[:,j] = WxNeumann[:,j] +
            wNeumannLeft[j-1] * WxNeumann[:,1]
        WxNeumann[:,end-(j-1)] = WxNeumann[:,end-(j-1)] +
            wNeumannRight[j-1] * WxNeumann[:,end]
        WhvDirichlet[:,j] = WhvDirichlet[:,j] +
            wDirichletLeft[j-1] * WhvDirichlet[:,1]
        WhvDirichlet[:,end-(j-1)] = WhvDirichlet[:,end-(j-1)] +
            wDirichletRight[j-1] * WhvDirichlet[:,end]
        WhvNeumann[:,j] = WhvNeumann[:,j] +
            wNeumannLeft[j-1] * WhvNeumann[:,1]
        WhvNeumann[:,end-(j-1)] = WhvNeumann[:,end-(j-1)] +
            wNeumannRight[j-1] * WhvNeumann[:,end]
    end

    WxDirichlet = WxDirichlet[:,2:end-1]
    WxNeumann = WxNeumann[:,2:end-1]
    WhvPeriodic = alp*dx^(phs-2) * WhvPeriodic
    WhvDirichlet = alp*dx^(phs-2) * WhvDirichlet[:,2:end-1]
    WhvNeumann = alp*dx^(phs-2) * WhvNeumann[:,2:end-1]
    
    #######################################################################
    # #TEST WxDirichlet and WxNeumann:
    # figure(1)
    # clf()
    # k = 3
    # plot(xi, WxDirichlet*sin.(k*pi*xi) - k*pi*cos.(k*pi*xi))
    # plot(xi, WxNeumann*cos.(k*pi*xi) + k*pi*sin.(k*pi*xi))
    # # plot(xi, fp(xi))
    # show()
    #######################################################################

    #Get block matrix and possibly plot eigenvalues:

    #Fast way:
    
    # A = vcat(hcat(WhvNeumann, WxDirichlet),
    #     hcat(WxNeumann, WhvDirichlet))

    # Aperiodic = vcat(hcat(WhvPeriodic, WxPeriodic),
    #     hcat(WxPeriodic, WhvPeriodic))

    #Other way:

    A = SparseMatrixCSC(zeros(2*ni, 2*ni))
    Aperiodic = SparseMatrixCSC(zeros(2*ni, 2*ni))

    Aperiodic[1:ni, ni+1:end] = WxPeriodic
    Aperiodic[ni+1:end, 1:ni] = WxPeriodic
    if plotEigenvalues
        eiPeriodicNohv = eigen(dt*Matrix(Aperiodic)).values
    end
    Aperiodic[1:ni,1:ni] = WhvPeriodic
    Aperiodic[ni+1:end, ni+1:end] = WhvPeriodic
    if plotEigenvalues
        eiPeriodic = eigen(dt*Matrix(Aperiodic)).values
    end

    A[1:ni, ni+1:end] = WxDirichlet
    A[ni+1:end, 1:ni] = WxNeumann
    if plotEigenvalues
        eiNohv = eigen(Matrix(dt*A)).values
    end
    A[1:ni, 1:ni] = WhvNeumann
    A[ni+1:end, ni+1:end] = WhvDirichlet
    if plotEigenvalues
        ei = eigen(Matrix(dt*A)).values
    end

    if plotEigenvalues
        a = min(minimum(real(eiNohv)), minimum(real(eiPeriodicNohv)))
        a = 1.05 * a
        b = -a
        c = -1.25
        d = 1.25
        figure(2)
        clf()
        subplot(131)
        plot(real(eiPeriodicNohv), imag(eiPeriodicNohv),
            marker=".", color="black", linestyle="none")
        title(@sprintf("periodic no HV, maxReal=%.2e",
            maximum(real(eiPeriodicNohv))))
        # axis(:image)
        axis([a,b,c,d])
        grid(:true)
        subplot(132)
        plot(real(eiNohv), imag(eiNohv), marker=".",
            color="black", linestyle="none")
        title(@sprintf("nonperiodic no HV, maxReal=%.2e",
            maximum(real(eiNohv))))
        # axis(:image)
        axis([a,b,c,d])
        grid(:true)
        a = minimum(real(ei))
        b = -0.05 * a
        a = 1.05 * a
        subplot(133)
        plot(real(ei), imag(ei),
            marker=".", color="black", linestyle="none")
        title(@sprintf("nonperiodic with HV, maxReal=%.2e",
            maximum(real(ei))))
        # axis(:image)
        axis([a,b,c,d])
        grid(:true)
        # subplot(144)
        # plot(real(eiPeriodic), imag(eiPeriodic),
        #     marker=".", color="black", linestyle="none")
        # title(@sprintf("p HV, maxReal=%.2e",
        #     maximum(real(eiPeriodic))))
        # axis(:image)
        # axis([a,b,c,d])
        # grid(:true)
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
            error(@sprintf("Unstable in time.  t=%g",t))
        end
        if plotSolution
            if mod(i, Int64(1/16/dx)) == 0
                figure(3)
                clf()
                subplot(121)
                plot(xi, U[1:ni])
                axis(:image)
                axis([-1,1,-1,1])
                title(@sprintf("rho, t=%1.2f",t))
                subplot(122)
                plot(xi, U[ni+1:end])
                axis(:image)
                axis([-1,1,-1,1])
                title(@sprintf("u, t=%1.2f",t))
                show()
                pause(.01)
            end
        end
    end

    if plotError
        figure(4)
        clf()
        subplot(121)
        plot(xi, U[1:ni] - U0[1:ni])
        title(@sprintf("error: rho, tf=%g",tf))
        subplot(122)
        plot(xi, U[ni+1:end])
        title(@sprintf("error: u, tf=%g",tf))
        show()
    end

end
