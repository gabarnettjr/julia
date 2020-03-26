########################################

using Plots

########################################

include("../gab/rk.jl")

########################################

#Initial time:
t0 = 0

#Final time:
tf = 5

#Number of time-steps:
N = 10

#RHS ODE function:
function f( t, y )
    return -y
end

#Initial condition:
y0 = 1

#Boolean:
exactSolutionKnown = true

#Exact solution:
function Y(t)
    return exp.(-t)
end

#Number of Runge-Kutta stages:
rkStages = 2

########################################

t = range( t0, stop = tf, length = N+1 )

y = zeros( size(t) )

y[1] = y0

h = ( tf - t0 ) / N

########################################

if rkStages == 1
    rk = rk1
elseif rkStages == 2
    rk = rk2
elseif rkStages == 3
    rk = rk3
elseif rkStages == 4
    rk = rk4
end

########################################

for n in 1 : N
    y[n+1] = rk( t[n], y[n], f, h )
end

########################################

p1 = plot( t, [y, Y(t)],
    title = "Numerical and Exact Solutions",
    label = ["numerical" "exact"],
    xlabel = "t", ylabel = "y" )

p2 = plot( t, y - Y(t), legend = false,
    title = "Difference" )

plot( p1, p2, layout = (1,2) )

########################################
