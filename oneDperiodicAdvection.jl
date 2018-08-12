using PyPlot
using Printf

include("rk.jl")
include("phs1.jl")

###########################################################################

pol = 5                             #highest degree polynomial in the basis
dx  = 1/128                                  #average spacing between nodes
dt  = dx/32                                               #constant delta t
tf  = 2                                                         #final time
ptb = 35                                      #node perturbation percentage
frq = 2                                     #frequency of initial sine wave
rks = 4                                                 #Runge-Kutta stages

###########################################################################

phs = 2*pol + 1
stc = 2*pol + 1

x = -1+dx/2 : dx : 1-dx/2
x = x + (-1 .+ 2*rand(Float64, length(x))) * dx * ptb/100

U0 = sin.(frq*pi*x)
U = copy(U0)

t = 0
nTimesteps = tf / dt

###########################################################################

#Get differentiation matrix with hyperviscosity:

Wx = getPeriodicDM(z=x, x=x, m=1,
    phs=phs, pol=pol, stc=stc, period=2)

Whv = getPeriodicDM(z=x, x=x, m=pol+1,
    phs=phs, pol=pol, stc=stc, period=2)

if pol == 1
    alp = 2.0 ^ -4
elseif pol == 3
    alp =  -2.0 ^ -7
elseif pol == 5
    alp =  2.0 ^ -10
elseif pol == 7
    alp = -2.0 ^ -13
else
    error("Choose odd integer for pol.")
end

W = -Wx + alp*dx^(pol-1) * Whv

###########################################################################

dUdt = zeros(size(x))
q1 = dUdt
q2 = zeros(size(x))
q3 = zeros(size(x))
q4 = zeros(size(x))

function odefun(t, U, dUdt)
    dUdt = W * U
    return dUdt
end

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
    if mod(i, Int64(1/4/dx)) == 0
        figure(1)
        plot(x, U)
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

figure(2)
plot(x, U-U0)
plot(x, zeros(size(x)), marker=".", linestyle="none")
title("error at tf, and nodes")
show()

###########################################################################
