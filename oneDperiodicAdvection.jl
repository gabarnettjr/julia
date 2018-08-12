using PyPlot
using Printf

include("rk.jl")
include("phsEigenvalues1d.jl")

###########################################################################

pol = 5                             #highest degree polynomial in the basis
dx  = 1/64                                   #average spacing between nodes
tf  = 2                                                         #final time
ptb = 35                                      #node perturbation percentage
frq = 5                                     #frequency of initial sine wave
rks = 4                                                 #Runge-Kutta stages
dt  = dx/8                                                #constant delta t

###########################################################################

e1, e2, x, W = PHSorFDeigs(ptb=ptb, frq=frq, dx=dx, dt=dt, pol=pol,
    showPlots=true)[1:4]

U0 = sin.(frq*pi*x)
U = copy(U0)

t = 0
nTimesteps = tf / dt

###########################################################################

dUdt = zeros(size(x))

function odefun(t, U, dUdt)
    dUdt = W * U
    return dUdt
end

###########################################################################

q1 = dUdt
q2 = zeros(size(x))
q3 = zeros(size(x))
q4 = zeros(size(x))

if rks == 3
    function rk!(t, U, odefun, dt)
        t, U = rk3!(t, U, odefun, dt, q1, q2)
        return t, U
    end
elseif rks == 4
    function rk!(t, U, odefun, dt)
        t, U = rk4!(t, U, odefun, dt, q1, q2, q3, q4)
        return t, U
    end
else
    error("Choose rks=3 or rks=4 for this problem.")
end

###########################################################################

for i in 1 : nTimesteps
    global t, U
    t, U = rk!(t, U, odefun, dt)
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

figure(3)
plot(x, U-U0)
plot(x, zeros(size(x)), marker=".", linestyle="none")
title("error at tf, and nodes")
show()

###########################################################################
