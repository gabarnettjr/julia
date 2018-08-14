using PyPlot
using Printf

include("rk.jl")
include("phsEigenvalues1d.jl")

###########################################################################

function timestepAndGetError(; showPlots=false, tf=2, ptb=35, frq=5,
    rks=4, pol=5, dx=1/32)
    #=
    Solve the simple advection equation u_t + u_x = 0 on periodic [-1,1]
    using either regular finite differences or polyharmonic spline finite
    differences.  If you choose an even number for pol, then it uses
    regular finite differences, but if you choose an odd number for pol,
    it uses a larger stencil size and include polyharmonic spline radial
    basis functions in the basis.
    =#
    # pol = 5                       #highest degree polynomial in the basis
    # dx  = 1/64                             #average spacing between nodes
    # tf  = 2                                                   #final time
    # ptb = 35                                #node perturbation percentage
    # frq = 5                               #frequency of initial sine wave
    # rks = 4                                           #Runge-Kutta stages
    dt = dx/8                                             #constant delta t
    
    x, W = PHSorFDeigs(ptb=ptb, frq=frq, dx=dx, dt=dt, pol=pol,
        showPlots=showPlots)
    
    # U0 = zeros(size(x))
    # ind = (x .< 0.5) .& (x .> -0.5)
    # U0[ind] = ones(size(x[ind]))
    U0 = cos.(frq*pi*x)
    U = copy(U0)
    
    t = 0
    nTimesteps = tf / dt
    
    dUdt = zeros(size(x))
    function odefun(t, U, dUdt)
        dUdt = W * U
        return dUdt
    end
    
    q1 = dUdt
    q2 = zeros(size(x))
    if rks == 4
        q3 = zeros(size(x))
        q4 = zeros(size(x))
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
        # if mod(i, Int64(1/4/dx)) == 0
        #     figure(2)
        #     plot(x, U)
        #     axis(:image)
        #     axis([-1,1,-1,1])
        #     title(@sprintf("numerical solution, t=%1.2f",t))
        #     show()
        #     pause(.01)
        #     if i != nTimesteps
        #         clf()
        #     end
        # end
    end
    
    if showPlots
        figure(3)
        plot(x, U-U0)
        # plot(x, U0)
        plot(x, zeros(size(x)), marker=".", linestyle="none")
        title(@sprintf("error at t=%g, and nodes", tf))
        show()
    end

    return maximum(abs.(U-U0))

end

###########################################################################
#TESTING:
tf  = 2
ptb = 40
frq = 5
rks = 4
pol = 5
#TEST 1 (eigenvalue plots):
timestepAndGetError(showPlots=true, tf=tf, ptb=ptb, frq=frq, rks=rks,
    pol=pol, dx=1/64)
#TEST 2 (convergence rate):
n = 6
maxErr = zeros(n)
for i in 4:4+n-1
    maxErr[i-3] = timestepAndGetError(showPlots=false, tf=tf, ptb=ptb,
        frq=frq, rks=rks, pol=pol, dx=1/2^i)
end
println(maxErr)
print(maxErr[1:end-1] ./ maxErr[2:end])
###########################################################################
