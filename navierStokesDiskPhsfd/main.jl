
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/disk/nodes.jl")

include("../packages/disk/navierStokes.jl")

include("../packages/phs2.jl")

include("../packages/rk.jl")

###############################################################################

# USER INPUT

# Number of layers of radial nodes on the unit disk (odd number)
layers = 17

# Set how much the nodes will be perturbed
ptb = 0

# Delta t
dt = 2^-10 / (layers - 1)

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

# Final time
tf = 1

###############################################################################

# Remove old files and remake the results directory
rm("results", recursive = true)
mkdir("results")

# Array of all times
t = range(0, stop=tf, step=dt)

# Get the nodes on the unit disk
x, y = makeRadialNodes(layers)

# Set the hyperviscosity coefficient
if K == 2
    a = 0
    # a = -1/2 / (layers - 1) ^ (2*K-1)
else
    error("K should be 2.")
end

# Index of the boundary nodes
bb = abs.(x .^ 2 .+ y .^ 2) .> (1 - 1e-6)

# Index of the interior nodes
ii = abs.(x .^ 2 .+ y .^ 2) .< (1 - 1e-6)

# Perturb the nodes by percentage ptb of node spacing
x, y = perturbNodes!(x, y, layers, ptb, bb)

# Expand the circle so it has a radius of 10 meters instead of 1
x = 10 * x
y = 10 * y

# Number of nodes total
n = length(x)

# Number of interior nodes
ni = length(x[ii])

# Number of boundary nodes
nb = length(x[bb])

# Get all of the DMs that will be needed for the ODE function
Wx = getDM(hcat(x,y)', hcat(x,y)', [1 0], phs, pol, stc, 0)
Wy = getDM(hcat(x,y)', hcat(x,y)', [0 1], phs, pol, stc, 0)
aWhv = a * getDM(hcat(x,y)', hcat(x,y)', [-1 -1], phs, pol, stc, K)

###############################################################################

# Initialize main solution array U

mu, k, lam, Cv, R = getConstants()

rho_0, u_0, v_0, e_0, p_0 = getInitialConditions(x, Cv, R)

U = zeros(length(x), 4)
U[:,1] = rho_0
U[:,2] = u_0
U[:,3] = v_0
U[:,4] = e_0

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
    ODEfunction!(t, U, dUdt, Wx, Wy, aWhv, rho_0, e_0, p_0, bb,
    mu, k, lam, Cv, R)

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

end

###############################################################################

# The main time-stepping loop

function mainLoop(t, U, odefun!, rk!, dt, tf, layers)

    frame = 0

    for i in 1 : Int(tf/dt) + 1

        if mod(i-1, 256) == 0
            # Print some info and save some things
            printAndSave(i, U, frame)
            frame = frame + 1
            # Time the Runge-Kutta update
            @time begin
                U = rk!(t[i], U, odefun!)
            end
        else
            # Just update the array U to the next time level
            U = rk!(t[i], U, odefun!)
        end
        
        # Stop running if the numerical solution blows up
        if (maximum(abs.(U[:,1])) > 10 || maximum(abs.(U[:,2])) > 50 ||
            maximum(abs.(U[:,3])) > 50 || maximum(abs.(U[:,4])) > 5*10^5)
            println("It blew up.")
            printAndSave(i, U, -999)
            break
        end

    end
    
end

mainLoop(t, U, odefun!, rk!, dt, tf, layers)

###############################################################################

# Save things that will be needed in the python plotting routine.

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

###############################################################################











