
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/nodes/disk/makeRadialNodes.jl")
include("../packages/nodes/disk/appendGhostNodes.jl")
include("../packages/nodes/connectedRectangles/zerosAndOnes.jl")
include("../packages/nodes/connectedDisks/makeNodes.jl")

include("../packages/eulerEquations/getConstants.jl")

include("../packages/eulerEquations/eulerian/freeSlipNoFlux.jl")
include("../packages/eulerEquations/eulerian/ODEfunction.jl")

include(string("../packages/eulerEquations/eulerian/connectedDisks/",
               "getInitialConditions.jl"))
include(string("../packages/eulerEquations/eulerian/connectedDisks/",
               "getIndices.jl"))
include(string("../packages/eulerEquations/eulerian/connectedDisks/",
               "tangentsNormals.jl"))

include("../packages/phs2.jl")

include("../packages/timeStepping/explicit/rk/rk3.jl")
include("../packages/timeStepping/explicit/ab/ab3.jl")

###############################################################################

# USER INPUT

# Flag to determine if the full code runs or if time-stepping is skipped
doTheTimeStepping = true

# Final time
tf = 2^-3

# The radius of the circle, in meters
radius = 10

# Number of layers of radial nodes on the unit disk (odd number)
layers = 33

# Delta t
dt = 2^-10 / (layers - 1)

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

# Get the nodes on the full domain (three connected disks)
x, y, dr, radius = makeNodes!(layers, radius)

# Set the hyperviscosity coefficient
if K == 2
    # a = 0
    a = -2^4 * dr ^ (2*K-1)
else
    error("K should be 2.")
end

###############################################################################

# Save x, y, and t for the python plotting routine.

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

io = open("./results/radius.txt", "w")
writedlm(io, radius, ' ')
close(io)

io = open("./results/dr.txt", "w")
writedlm(io, dr, ' ')
close(io)

io = open("./results/t.txt", "w")
writedlm(io, t, ' ')
close(io)

###############################################################################

# Get and save all of the indices for the various layers of nodes.

nb = 6 * (layers - 3)
dth = 2 * pi / nb
indA, indB, indC = getIndices(x, y, dr, radius, 4 * dth + 1e-3)

Tx, Ty, Nx, Ny = tangentsNormals((x[indB].+x[indC])./2,
                                 (y[indB].+y[indC])./2, radius, dr)

io = open("./results/indA.txt", "w")
writedlm(io, indA, ' ')
close(io)

io = open("./results/indB.txt", "w")
writedlm(io, indB, ' ')
close(io)

io = open("./results/indC.txt", "w")
writedlm(io, indC, ' ')
close(io)

io = open("./results/doTheTimeStepping.txt", "w")
writedlm(io, doTheTimeStepping, ' ')
close(io)

###############################################################################

# Get all of the DMs that will be needed for the ODE function

Wx, tr, ind_nn =
     getDM(hcat(x,y)', hcat(x,y)', [1 0], phs, pol, stc, 0)[[1,3,4]]

Wy = getDM(hcat(x,y)', hcat(x,y)', [0 1], phs, pol, stc, 0;
           tree = tr, idx = ind_nn)[1]

aWhv = a * getDM(hcat(x,y)', hcat(x,y)', [-1 -1], phs, pol, stc, K;
                 tree = tr, idx = ind_nn)[1]

###############################################################################

# Initialize main solution array U and save initial conditions.

Cv, R = getConstants()

rho_0, u_0, v_0, e_0, p_0 = getInitialConditions(x, Cv, R, radius, dr)

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

U = zeros(length(x), 4)
U[:,1] = rho_0
U[:,2] = u_0
U[:,3] = v_0
U[:,4] = e_0

###############################################################################

# Initialize dummy arrays to be used in Runge-Kutta and Adams-Bashforth

q1 = zeros(size(U))                                             #needed for RK3
q2 = zeros(size(U))                                             #needed for RK3

f1 = zeros(size(U))                                             #needed for AB3
f2 = zeros(size(U))                                             #needed for AB3
f3 = zeros(size(U))                                             #needed for AB3

###############################################################################

# Define an ODE function which can be passed to rk

odefun! = (t, U, dUdt) -> ODEfunction!(t, U, dUdt,
                                       Wx, Wy, aWhv, rho_0, e_0,
                                       Tx, Ty, Nx, Ny, indA, indB, indC,
                                       Cv, R)

###############################################################################

function printAndSave(i, U, frame)

    @printf("frame = %.0f,  t = %.5f\n", frame, t[i])
    @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f,  maxE = %.5f\n", 
            maximum(U[:,1]-rho_0), maximum(U[:,2]), maximum(U[:,3]),
            maximum(U[:,4]-e_0))
    @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f,  minE = %.5f\n\n", 
            minimum(U[:,1]-rho_0), minimum(U[:,2]), minimum(U[:,3]),
            minimum(U[:,4]-e_0))

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

function mainLoop(U, f1, f2, f3)

    # Save initial stuff
    frame = 0
    U = freeSlipNoFlux!(U, Tx, Ty, Nx, Ny, indA, indB, indC)
    stop = printAndSave(1, U, frame)
    frame = 1

    # Do two steps of Runge-Kutta, to get first three time-levels
    f1 = odefun!(t[1], U, f1)
    U = rk3!(t[1], U, odefun!, dt, q1, q2)
    f2 = odefun!(t[2], U, f2)
    U = rk3!(t[2], U, odefun!, dt, q1, q2)
    f3 = odefun!(t[3], U, f3)

    for i in 3 : Int(tf/dt) + 1

        if mod(i-1, 2^7) == 0
            U = freeSlipNoFlux!(U, Tx, Ty, Nx, Ny, indA, indB, indC)
            # Print some info and save some things
            stop = printAndSave(i, U, frame)
            if stop
                break
            end
            frame = frame + 1
            # Time the update
            @time begin
                U, f1, f2, f3 = ab3!(t[i], U, odefun!, dt, f1, f2, f3)
            end
        else
            # Just step forward in time
            U, f1, f2, f3 = ab3!(t[i], U, odefun!, dt, f1, f2, f3)
        end

    end
    
end

if doTheTimeStepping
    mainLoop(U, f1, f2, f3)
end

###############################################################################











