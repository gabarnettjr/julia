
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/nodes/connectedSquares/zerosAndOnes.jl")
include("../packages/nodes/connectedSquares/adjacency.jl")
include("../packages/nodes/connectedSquares/refine.jl")
include("../packages/nodes/connectedSquares/makeNodes.jl")
include("../packages/nodes/connectedSquares/getBooleans.jl")
include("../packages/nodes/connectedSquares/getIndices.jl")

include("../packages/eulerEquations/getConstants.jl")

include(string(("../packages/eulerEquations/eulerian/connectedSquares/",
                "matsToVecs.jl"))
include(string(("../packages/eulerEquations/eulerian/connectedSquares/",
                "vecsToMats.jl"))
# include(string("../packages/eulerEquations/eulerian/connectedSquares/",
#                "getInitialConditions.jl"))
include(string("../packages/eulerEquations/eulerian/connectedSquares/",
               "freeSlipNoFlux.jl"))
# include(string("../packages/eulerEquations/eulerian/connectedSquares/",
#                "ODEfunction.jl"))

include("../packages/phs2.jl")

include("../packages/timeStepping/explicit/rk/rk3.jl")
include("../packages/timeStepping/explicit/ab/ab3.jl")

###############################################################################

# USER INPUT

# Final time
tf = 2^-3

# Choose which domain to use
domain = 1

# The level of refinement (lowest is 1, dx = 1 / refinement)
refinement = 1

# Delta t
dt = 2^-10 / refinement

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

# Get the matrix of zeros and ones that defines the room shape
z1 = zerosAndOnes(domain = 1)

# Get the matrix of numbers that defines the adjacency of each square
A = adjacency(z1)

# Refine to the requested refinement level
z1, A = refine!(z1, A, refinement)

# Define the x and y coordinates (x and y are matrices, meshgrid style)
x, y, dx = makeNodes(z1, refinement)

# Get the BitArray matrices to extract certain nodes
bool_all, bool_noGhost = getBooleans(A)

# Get all the indices that might be used to enforce the boundary conditions
i_n1, j_n1, i_n2, j_n2, i_n3, j_n3, i_n4, j_n4,
i_n12, j_n12, i_n13, j_n13, i_n14, j_n14, i_n23, j_n23,
i_n24, j_n24, i_n34, j_n34,
i_n123, j_n123, i_n124, j_n124, i_n134, j_n134, i_n234, j_n234,
i_n1234, j_n1234 = getIndices(A)

# Set the hyperviscosity coefficient
if K == 2
    # a = 0
    a = -2^4 * dx ^ (2*K-1)
else
    error("K should be 2.")
end

###############################################################################

# Get and save all of the indices for the various types of nodes.

io = open("./results/indA.txt", "w")
writedlm(io, indA, ' ')
close(io)

io = open("./results/indB.txt", "w")
writedlm(io, indB, ' ')
close(io)

io = open("./results/indC.txt", "w")
writedlm(io, indC, ' ')
close(io)

###############################################################################

# Get all of the DMs that will be needed for the ODE function

Wx, tr, ind_nn = getDM(hcat(x[bool_noGhost], y[bool_noGhost])',
                       hcat(x[bool_all], y[bool_all])',
                       [1 0], phs, pol, stc, 0)[[1,3,4]]

Wy = getDM(hcat(x[bool_noGhost], y[bool_noGhost])',
           hcat(x[bool_all], y[bool_all])',
           [0 1], phs, pol, stc, 0;
           tree = tr, idx = ind_nn)[1]

aWhv = a * getDM(hcat(x[ind_noGhost], y[ind_noGhost])', 
                 hcat(x[ind_all], y[ind_all])',
                 [-1 -1], phs, pol, stc, K;
                 tree = tr, idx = ind_nn)[1]

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

mainLoop(U, f1, f2, f3)

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











