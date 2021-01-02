
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/nodes/tangents.jl")

include("../packages/eulerEquations/getConstants.jl")

###############################################################################

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "getInitialConditions.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "getNodes.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "getIndices.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "removeNormalComponent.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "fanPrime.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "ODEfunction.jl"))

include(string("../packages/eulerEquations/eulerian/boundaryDisk/",
               "simpleFixer.jl"))

###############################################################################

include("../packages/phs2.jl")

# include("../packages/timeStepping/explicit/rk/rk3.jl")
include("../packages/timeStepping/explicit/rk/rk4.jl")
# include("../packages/timeStepping/explicit/ab/ab3.jl")

###############################################################################

# USER INPUT

# Final time (power of 2)
tf = 2^5

# The radius of the circle in meters (power of 2)
radius = 2^3

# Number of layers of radial layers on the disk (power of 2 PLUS ONE)
layers = 2^5 + 1

# The number of times to save (power of 2 PLUS ONE)
numSaves = 2^10 + 1

# Approximate node spacing
dr = radius / (layers - 1)

# Delta t (power of 2 times dr)
dt = 2^-10 * dr

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

# Array of times when things are saved
tSaves = range(0, stop=tf, length = numSaves)

# Get the nodes on the disk, with the wall built in
x, y, xNon, yNon, s, doorWidth, outAngle, dth =
    getNodes(layers, radius)

# Get the index of the different types of nodes
bb, ff, ii = getIndices(x, y, dr, s, doorWidth, outAngle, dth)

# Set the hyperviscosity coefficient
if K == 2
    # a = 0
    a = -2^2 * dr ^ (2*K-1)
else
    error("K should be 2.")
end

# Get the unit tangent vector at each boundary node
Tx, Ty = tangents(x[bb], y[bb])

# Physical constants having to do with the fluid (dry air)
Cv, R = getConstants()

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

io = open("./results/dr.txt", "w")
writedlm(io, dr, ' ')
close(io)

io = open("./results/Cv.txt", "w")
writedlm(io, Cv, ' ')
close(io)

io = open("./results/R.txt", "w")
writedlm(io, R, ' ')
close(io)

io = open("./results/tSaves.txt", "w")
writedlm(io, tSaves, ' ')
close(io)

io = open("./results/bb.txt", "w")
writedlm(io, bb, ' ')
close(io)

io = open("./results/ff.txt", "w")
writedlm(io, ff, ' ')
close(io)

io = open("./results/xNon.txt", "w")
writedlm(io, xNon, ' ')
close(io)

io = open("./results/yNon.txt", "w")
writedlm(io, yNon, ' ')
close(io)

###############################################################################

# Get all of the DMs that will be needed in the ODE function

Wx, tr, ind_nn = getDM(hcat(x,y)', hcat(x,y)',
                       [1 0], phs, pol, stc, 0)[[1,4,5]]

Wy = getDM(hcat(x,y)', hcat(x,y)',
           [0 1], phs, pol, stc, 0;
           tree = tr, idx = ind_nn)[1]

aWhv = a .* getDM(hcat(x,y)', hcat(x,y)',
                  [-1 -1], phs, pol, stc, K;
                  tree = tr, idx = ind_nn)[1]

###############################################################################

# Initialize main solution array U and save initial conditions.

rho_0, e_0, q_0, p_0 = getInitialConditions(x, Cv, R, s, radius, dr)

# Get initial averages for density, tracer density, and energy
rho_0_avg = sum(rho_0) / length(rho_0)
rhoq_0_avg = sum(rho_0 .* q_0) / length(q_0)
e_0_avg = sum(e_0) / length(e_0)

io = open("./results/rho_0.txt", "w")
writedlm(io, rho_0, ' ')
close(io)

io = open("./results/e_0.txt", "w")
writedlm(io, e_0, ' ')
close(io)

io = open("./results/q_0.txt", "w")
writedlm(io, q_0, ' ')
close(io)

U = zeros(length(x), 4)
U[:,1] = rho_0
U[:,2] .= 0
U[:,3] .= 0
U[:,4] = e_0
q = copy(q_0)

###############################################################################

# Initialize dummy arrays to be used in Runge-Kutta and Adams-Bashforth

q1 = zeros(size(U))                                             #needed for RK3
q2 = zeros(size(U))                                             #needed for RK3
q3 = zeros(size(U))                                             #needed for RK4
q4 = zeros(size(U))                                             #needed for RK4

f1 = zeros(size(U))                                             #needed for AB3
f2 = zeros(size(U))                                             #needed for AB3
f3 = zeros(size(U))                                             #needed for AB3

###############################################################################

# Define an ODE function which can be passed to rk

odefun! = (t, U, dUdt) -> ODEfunction!(t, U, dUdt,
                                       Wx, Wy, aWhv, rho_0, e_0, p_0,
                                       Tx, Ty, bb, ff,
                                       Cv, R)

###############################################################################

function printAndSave!(i, U, q, xt, yt, frame,
                       maxVel, minRho, maxRho, minU, maxU, minV, maxV,
                       minE, maxE, minQ, maxQ, minP, maxP)
    
    maxVel = vcat(maxVel, maximum(sqrt.(U[:,2] .^ 2 .+ U[:,3] .^ 2)))
    minRho = vcat(minRho, minimum(U[:,1] - rho_0))
    maxRho = vcat(maxRho, maximum(U[:,1] - rho_0))
    minU   = vcat(minU,   minimum(U[:,2]))
    maxU   = vcat(maxU,   maximum(U[:,2]))
    minV   = vcat(minV,   minimum(U[:,3]))
    maxV   = vcat(maxV,   maximum(U[:,3]))
    minE   = vcat(minE,   minimum(U[:,4] - e_0))
    maxE   = vcat(maxE,   maximum(U[:,4] - e_0))
    minQ   = vcat(minQ,   minimum(q))
    maxQ   = vcat(maxQ,   maximum(q))
    minP   = vcat(minP,   minimum(U[:,1] .* R .* (U[:,4] ./ Cv) - p_0))
    maxP   = vcat(maxP,   maximum(U[:,1] .* R .* (U[:,4] ./ Cv) - p_0))

    @printf("frame = %.0f,  t = %.5f\n", frame, t[i])
    @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f,  maxE = %.5f,  maxQ = %.5f\n", 
            maxRho[end],
            maxU[end],
            maxV[end],
            maxE[end],
            maxQ[end])
    @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f,  minE = %.5f,  minQ = %.5f\n\n", 
            minRho[end],
            minU[end],
            minV[end],
            minE[end],
            minQ[end])

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

    io = open(@sprintf("./results/q_%04d.txt",frame), "w")
    writedlm(io, q, ' ')
    close(io)

    io = open(@sprintf("./results/xt_%04d.txt",frame), "w")
    writedlm(io, xt, ' ')
    close(io)

    io = open(@sprintf("./results/yt_%04d.txt",frame), "w")
    writedlm(io, yt, ' ')
    close(io)

    stop = false

    # Stop running if the numerical solution blows up
    if (maximum(abs.(U[:,1])) > 10 ||
        maximum(abs.(U[:,2])) > 50 ||
        maximum(abs.(U[:,3])) > 50 ||
        maximum(abs.(U[:,4])) > 5*10^5 ||
        !iszero(isnan.(U)))
        println("It blew up.")
        stop = true
    end

    return stop, maxVel, minRho, maxRho, minU, maxU, minV, maxV, minE, maxE,
           minQ, maxQ, minP, maxP

end

###############################################################################

# The main time-stepping loop

function mainLoop!(U, q, f1, f2, f3)

    maxVel = []

    minRho = [];  maxRho = []
    minU   = [];  maxU   = []
    minV   = [];  maxV   = []
    minE   = [];  maxE   = []
    minQ   = [];  maxQ   = []
    minP   = [];  maxP   = []

    # Save initial stuff
    frame = 0
    stop, maxVel, minRho, maxRho, minU, maxU, minV, maxV, minE, maxE,
    minQ, maxQ, minP, maxP = printAndSave!(1, U, q, x, y, frame, maxVel,
                               minRho, maxRho, minU, maxU, minV, maxV,
                               minE, maxE, minQ, maxQ, minP, maxP)
    frame = 1

    # Testing out some methods for keeping track of moving points.  (xt,yt)
    # are the locations of moving particles, which will follow the flow.
    DT = tf / (numSaves - 1)
    xt = copy(x)
    yt = copy(y)
    u0 = copy(U[:,2])
    v0 = copy(U[:,3])

    # Nodes that will move for only one large time step for SL transport
    xq = copy(x)
    yq = copy(y)

    # Do two steps of Runge-Kutta for prog vars, to get first three time-levels
    f1 = odefun!(t[1], U, f1)
    U = rk4!(t[1], U, odefun!, dt, q1, q2, q3, q4)
    f2 = odefun!(t[2], U, f2)
    U = rk4!(t[2], U, odefun!, dt, q1, q2, q3, q4)
    f3 = odefun!(t[3], U, f3)

    for i in 3 : Int(tf/dt) + 1

        if mod(i-1, Int(tf/dt/(numSaves-1))) == 0
            # Let the moving points take a leap forward to the present time
            W = getDM(hcat(xt,yt)', hcat(x,y)',
                      [0 0], phs, pol, stc, 0;
                      tree = tr)[1]
            xHalf = xt .+ (DT/2) .* (W * u0)
            yHalf = yt .+ (DT/2) .* (W * v0)
            uHalf = (u0 .+ U[:,2]) ./ 2
            vHalf = (v0 .+ U[:,3]) ./ 2
            W = getDM(hcat(xHalf,yHalf)', hcat(x,y)',
                      [0 0], phs, pol, stc, 0;
                      tree = tr)[1]
            xt = xt .+ DT .* (W * uHalf)
            yt = yt .+ DT .* (W * vHalf)
            # Put wayward points back in the domain randomly
            r = sqrt.(xt .^ 2 .+ yt .^ 2)
            ind_out = r .> radius
            ran = rand(length(r[ind_out]))
            xt[ind_out] = xt[ind_out] ./ r[ind_out] .* radius .* ran
            yt[ind_out] = yt[ind_out] ./ r[ind_out] .* radius .* ran
            # Do one step of semi-lagrangian tracer transport
            xHalf = x .+ (DT/2) .* u0
            yHalf = y .+ (DT/2) .* v0
            uHalf = (u0 .+ U[:,2]) ./ 2
            vHalf = (v0 .+ U[:,3]) ./ 2
            W = getDM(hcat(xHalf,yHalf)', hcat(x,y)',
                      [0 0], phs, pol, stc, 0;
                      tree = tr, idx = ind_nn)[1]
            xq = x .+ DT .* (W * uHalf)
            yq = y .+ DT .* (W * vHalf)
            W = getDM(hcat(x,y)', hcat(xq,yq)',
                      [0 0], phs, pol, stc, 0)[1]
            q = W * q
            # Shift saved velocity to the new time level
            u0 = U[:,2]
            v0 = U[:,3]
            # Apply a fixer for mass, tracer mass, and total energy
            U, q = simpleFixer!(U, q, rho_0_avg, rhoq_0_avg, e_0_avg, ii)
            # Print some info and save some things
            stop, maxVel, minRho, maxRho, minU, maxU, minV, maxV, minE, maxE,
            minQ, maxQ, minP, maxP = printAndSave!(i, U, q, xt, yt, frame,
                                       maxVel, minRho, maxRho, minU, maxU,
                                       minV, maxV, minE, maxE, minQ, maxQ,
                                       minP, maxP)
            # Time the update
            @time begin
                # U, f1, f2, f3 = ab3!(t[i], U, odefun!, dt, f1, f2, f3)
                # U = rk3!(t[i], U, odefun!, dt, q1, q2)
                U = rk4!(t[i], U, odefun!, dt, q1, q2, q3, q4)
            end
            # Stop if it blew up
            if stop
                break
            end
            frame = frame + 1
        else
            # Just step forward in time
            # U, f1, f2, f3 = ab3!(t[i], U, odefun!, dt, f1, f2, f3)
            # U = rk3!(t[i], U, odefun!, dt, q1, q2)
            U = rk4!(t[i], U, odefun!, dt, q1, q2, q3, q4)
        end

    end
    
    io = open("./results/maxVel.txt", "w")
    writedlm(io, maximum(maxVel), ' ')
    close(io)

    io = open("./results/minRho.txt", "w")
    writedlm(io, minimum(minRho), ' ')
    close(io)

    io = open("./results/maxRho.txt", "w")
    writedlm(io, maximum(maxRho), ' ')
    close(io)

    io = open("./results/minU.txt", "w")
    writedlm(io, minimum(minU), ' ')
    close(io)

    io = open("./results/maxU.txt", "w")
    writedlm(io, maximum(maxU), ' ')
    close(io)

    io = open("./results/minV.txt", "w")
    writedlm(io, minimum(minV), ' ')
    close(io)

    io = open("./results/maxV.txt", "w")
    writedlm(io, maximum(maxV), ' ')
    close(io)

    io = open("./results/minE.txt", "w")
    writedlm(io, minimum(minE), ' ')
    close(io)

    io = open("./results/maxE.txt", "w")
    writedlm(io, maximum(maxE), ' ')
    close(io)

    io = open("./results/minQ.txt", "w")
    writedlm(io, minimum(minQ), ' ')
    close(io)

    io = open("./results/maxQ.txt", "w")
    writedlm(io, maximum(maxQ), ' ')
    close(io)

    io = open("./results/minP.txt", "w")
    writedlm(io, minimum(minP), ' ')
    close(io)

    io = open("./results/maxP.txt", "w")
    writedlm(io, maximum(maxP), ' ')
    close(io)

end

mainLoop!(U, q, f1, f2, f3)

###############################################################################













