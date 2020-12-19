
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/disk/nodes.jl")

include("../packages/disk/eulerEquationsSemiLagrangian.jl")

include("../packages/phs2.jl")

include("../packages/rk.jl")

###############################################################################

# USER INPUT

# Final time
tf = 1

# The radius of the circle, in meters
radius = 10

# Number of layers of radial nodes on the unit disk (odd number)
layers = 33

# Delta t
dt = 2^-5 / (layers - 1)

# Runge-Kutta stages
rkstages = 3

# Exponent in the polyharmonic spline function
phs = 5

# Highest degree of polynomial to include in the basis
pol = 2

# Stencil size
stc = 19

# Hyperviscosity exponent
K = 2

###############################################################################

# Remove old files and remake the results directory
rm("results", recursive = true)
mkdir("results")

# Array of all times
t = range(0, stop=tf, step=dt)

# Get the nodes on the unit disk
x, y = makeRadialNodes(layers - 2)
dr = 1 / (layers - 3)

# Modify the nodes a bit so there are matching layers near the boundary:
x, y = appendGhostNodes!(x, y, dr)

# Scale the nodes so the radius is correct (nearly):
x = radius .* x
y = radius .* y
dr = sqrt((x[1] - x[2])^2 + (y[1] - y[2])^2)
radius = radius + 3 ./ 2 .* dr

# Set the hyperviscosity coefficient
if K == 2
    a = 0
    # a = -2^-1 * dr ^ (2*K-1)
else
    error("K should be 2.")
end

###############################################################################

# Get and save all of the indices for the various layers of nodes.
ind_ghost, ind_noGhost, ind_outer, ind_fan = getIndices(x, y, dr, radius)

io = open("./results/ind_ghost.txt", "w")
writedlm(io, ind_ghost, ' ')
close(io)

io = open("./results/ind_noGhost.txt", "w")
writedlm(io, ind_noGhost, ' ')
close(io)

io = open("./results/ind_outer.txt", "w")
writedlm(io, ind_outer, ' ')
close(io)

io = open("./results/ind_fan.txt", "w")
writedlm(io, ind_fan, ' ')
close(io)

###############################################################################

# Initialize main solution array U and save initial conditions.

Cv, R = getConstants()

rho_0, u_0, v_0, e_0, p_0 = getInitialConditions(x, Cv, R)

io = open("./results/rho_0.txt", "w")
writedlm(io, rho_0, ' ')
close(io)

io = open("./results/u_0.txt", "w")
writedlm(io, u_0, ' ')
close(io)

io = open("./results/v_0.txt", "w")
writedlm(io, v_0, ' ')
close(io)

io = open("./results/e_0.txt", "w")
writedlm(io, e_0, ' ')
close(io)

U = zeros(length(x), 6)
U[:,1] = rho_0
U[:,2] = u_0
U[:,3] = v_0
U[:,4] = e_0
U[:,5] = x
U[:,6] = y

###############################################################################

# Initialize dummy arrays to be used in Runge-Kutta

q1 = zeros(size(U))                                        #needed in all cases
q2 = zeros(size(U))                                     #needed for rk3 and rk4
q3 = []
q4 = []
if rkstages == 4
    q3 = zeros(size(U))                                         #needed for rk4
    q4 = zeros(size(U))                                         #needed for rk4
end

###############################################################################

# Define an ODE function which can be passed to rk

odefun! = (t, U, dUdt) ->
    ODEfunction!(t, U, dUdt,
    radius, rho_0, e_0, p_0, phs, pol, stc, K,
    ind_ghost, ind_noGhost, ind_outer, ind_fan,
    Cv, R, a)
    
###############################################################################

# Define the Runge-Kutta function based on the number of stages

if rkstages == 3
	rk! = (t, U, odefun!) -> rk3!(t, U, odefun!, dt, q1, q2)
elseif rkstages == 4
	rk! = (t, U, odefun!) -> rk4!(t, U, odefun!, dt, q1, q2, q3, q4)
else
    error("rkstages should be 3 or 4 please.")
end

###############################################################################

function printAndSave(i, U, frame)

    @printf("i = %.0f,  t = %.5f\n", i, t[i])
    @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f,  maxE = %.5f\n", 
            maximum(U[:,1].-rho_0), maximum(U[:,2]), maximum(U[:,3]),
            maximum(U[:,4].-e_0))
    @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f,  minE = %.5f\n\n", 
            minimum(U[:,1].-rho_0), minimum(U[:,2]), minimum(U[:,3]),
            minimum(U[:,4].-e_0))

    io = open(@sprintf("./results/rho_%04d.txt",frame), "w")
    writedlm(io, U[:,1], ' ')
    close(io)

    io = open(@sprintf("./results/u_%04d.txt",frame), "w")
    writedlm(io, U[:,2], ' ')
    close(io)

    io = open(@sprintf("./results/v_%04d.txt",frame), "w")
    writedlm(io, U[:,3], ' ')
    close(io)
    
    io = open(@sprintf("./results/e_%04d.txt",frame), "w")
    writedlm(io, U[:,4], ' ')
    close(io)

    stop = false

    # Stop running if the numerical solution blows up
    if (maximum(abs.(U[:,1])) > 10 || maximum(abs.(U[:,2])) > 50 ||
        maximum(abs.(U[:,3])) > 50 || maximum(abs.(U[:,4])) > 5*10^5 ||
        !iszero(isnan.(U)))
        println("It blew up.")
        stop = true
    end

    return stop

end

###############################################################################

# The main time-stepping loop

function mainLoop(t, U, odefun!, rk!, dt, tf, layers, x, y,
ind_ghost, ind_noGhost, ind_outer, ind_fan,
phs, pol, stc, e_0)

    frame = 0

    for i in 1 : Int(tf/dt) + 1

        if mod(i-1, 2^0) == 0
            U = freeSlipNoFlux!(U, ind_ghost, ind_noGhost, ind_outer, radius,
                                phs, pol, stc, e_0)
            # U = fan!(t[i], U, ind_fan, phs, pol, stc)
            # Print some info and save some things
            stop = printAndSave(i, U, frame)
            if stop
                break
            end
            frame = frame + 1
            # Time the update
            @time begin
                U = rk!(t[i], U, odefun!)
            end
        else
            # Just step forward in time
            U = rk!(t[i], U, odefun!)
        end

        # Interpolate everything back to the original nodes:
        U = interp!(U, x, y, phs, pol, stc)

    end
    
end

mainLoop(t, U, odefun!, rk!, dt, tf, layers, x, y,
         ind_ghost, ind_noGhost, ind_outer, ind_fan,
         phs, pol, stc, e_0)

###############################################################################

# Save things that will be needed in the python plotting routine.

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

io = open("./results/radius.txt", "w")
writedlm(io, radius, ' ')
close(io)

io = open("./results/t.txt", "w")
writedlm(io, t, ' ')
close(io)

###############################################################################











