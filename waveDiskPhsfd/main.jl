#####################################################################

using Printf

#####################################################################

include("../packages/disk/nodes.jl")

include("../packages/disk/waveEquation.jl")

include("../packages/rk.jl")

#####################################################################

# USER INPUT

# Number of layers of radial nodes on the unit disk
layers = 33

# Delta t
dt = 1/2 * 1/(layers - 1)

# Runge-Kutta stages
rkstages = 3

# Exponent in the polyharmonic spline function
phs = 5

# Highest degree of polynomial to include in the basis
pol = 3

# Stencil size
stc = 37

# Hyperviscosity exponent
K = 2

# Wave speed
c = 1.

# Final time
tf = 10

#####################################################################

t = range(0, stop=tf, step=dt)

x, y = makeRadialNodes(layers)

if K == 2
    a = 0
    # a = -2^(-5) * c * 1 / (layers - 1) ^ (2*K-1)
else
    error("Still need to implement other cases.")
end

# Index of the interior nodes
ii = abs.(x .^ 2 .+ y .^ 2) .< (1 - 1e-6)

cWx, cWy, aWhv = getAllDMs(c, hcat(x,y)', ii, phs, pol, stc, K, a)

#####################################################################

# Initialize main solution array U
U = zeros(length(x), 3)
U = initialCondition!(x, y, U)

# Initialize dummy arrays to be used in Runge-Kutta
q1 = zeros(size(U))
if rkstages > 2
    q2 = zeros(size(U))
end
if rkstages > 3
    q3 = zeros(size(U))
    q4 = zeros(size(U))
end

#####################################################################

# Define an ODE function which can be passed to rk.jl

function odefun!(t, U, dUdt)

    global cWx, cWy, aWhv, ii

    dUdt = ODEfunction!(t, U, dUdt, cWx, cWy, aWhv, ii)

    return dUdt

end

#####################################################################

# Define the Runge-Kutta function based on the number of stages

if rkstages == 2

    function rk!(t, U, odefun!, dt)
        return rk2!(t, U, odefun!, dt, q1)
    end

elseif rkstages == 3

    function rk!(t, U, odefun!, dt)
        return rk3!(t, U, odefun!, dt, q1, q2)
    end

elseif rkstages == 4

    function rk!(t, U, odefun!, dt)
        return rk4!(t, U, odefun!, dt, q1, q2, q3, q4)
    end

else

    error("rkstages should be 2, 3, or 4 please.")

end

#####################################################################

# The main time-stepping loop

numSaves = 0

for i in 1 : Int(tf/dt)+1

    global rk!, odefun!, rkstages, U, t, dt, q1, q2, q3, q4, numSaves

    if mod(i-1, (layers-1)) == 0
        @printf("i = %.0f,  t = %.5f\n", i, t[i])
        @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f\n", 
                maximum(U[:,1]), maximum(U[:,2]), maximum(U[:,3]))
        @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f\n\n", 
                minimum(U[:,1]), minimum(U[:,2]), minimum(U[:,3]))
        numSaves = numSaves + 1
        io = open(@sprintf("./results/rho_%04d.txt",numSaves), "w")
        writedlm(io, U[:,1], ' ')
        close(io)
    end

    U = rk!(t[i], U, odefun!, dt)

    if maximum(abs.(U)) > 20
        println("It blew up.")
        break
    end

end

#####################################################################










