
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

include("../packages/nodes/getFloorplanMatrix.jl")
include("../packages/nodes/expandFloorplanMatrix.jl")
include("../packages/nodes/getFloorplanNodesAndIndices.jl")
include("../packages/nodes/findCornersAndSides.jl")
include("../packages/nodes/tangents.jl")
include("../packages/nodes/countIntegrationSquares.jl")

include("../packages/eulerEquations/getConstants.jl")

include("../packages/vectors/removeNormalComponent.jl")

include("../packages/phs2.jl")

include("../packages/timeStepping/explicit/rk/rk3.jl")

###############################################################################

include("getInitialConditions.jl")
include("simpleFixer.jl")
include("fanPrime.jl")
include("ODEfunction.jl")

###############################################################################

# USER INPUT

# Final time (power of 2)
tf = 2^4

# The type of floorplan
floorplan = "vinny5"

# Speed of air blown by the fan (m/s)
fanSpeed = 10

# Approximate node spacing in meters.  The wall width is always 1, but dr
# controls how many layers wide each wall is.  dr=1/7 => 6 layers wide.
# This should be one divided by an odd number (5, 7, 9, ...)
k = 5
dr = 1/k
overlap = Int((k-1)/2) * dr

# Delta t in seconds (power of 2 times dr)
dt = 2^-11 * dr

# Array of all times in seconds
t = range(0, stop=tf, step=dt)

# The number of times to save (power of 2 PLUS ONE)
numSaves = 2^(Int(log2(tf))+6) + 1

# Exponent in the polyharmonic spline function
phs = 5

# Highest degree of polynomial to include in the basis
pol = 2

# Stencil size
stc = 15

# Hyperviscosity exponent
K = 2

###############################################################################

# Remove old files and remake the results directory
rm("results", recursive = true)
mkdir("results")

# Array of times when things are saved
tSaves = range(0, stop=tf, length = numSaves)

# Get the nodes on the domain, with the door built in
M_0, len, wid = getFloorplanMatrix(plan = floorplan)
M = expandFloorplanMatrix(M_0, Int(1/dr))
x, y, xNon, yNon, bb, ff, ii, qq =
    getFloorplanNodesAndIndices(M, len, wid, dr)

# Find corner and side nodes to be used in mass fixers
cs1, cs2, cs3, cs4 = findCornersAndSides(x, y, dr)

# Count the number of small squares, to be used in the mass fixers
numSquares = countIntegrationSquares(x, y, len, wid, dr)

# Set the hyperviscosity coefficient
if K == 2
    # a = 0
    a = -2^3 * dr ^ (2*K-1)
elseif K == 3
    a = 2^0 * dr ^ (2*K-1)
else
    error("K should be 2 or 3.")
end

# Get the unit tangent vector at each boundary node
Tx, Ty = tangents(x[bb], y[bb])

# # Repeated version of Tx and Ty, to be used in freeSlipNoFlux!()
# TX = repeat(Tx', stc)
# TY = repeat(Ty', stc)

# Physical constants having to do with the fluid (dry air)
Cp, Cv, R = getConstants()

###############################################################################

# Save things that will be needed in the python plotting routine.

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

io = open("./results/dr.txt", "w")
writedlm(io, dr, ' ')
close(io)

io = open("./results/Cp.txt", "w")
writedlm(io, Cp, ' ')
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

Wx, wx, jj, tr, ind_nn = getDM(hcat(x,y)', hcat(x,y)',
                               [1 0], phs, pol, stc, 0)

Wy, wy = getDM(hcat(x,y)', hcat(x,y)',
               [0 1], phs, pol, stc, 0;
               tree = tr, idx = ind_nn)[[1,2]]

aWhv, awhv = a .* getDM(hcat(x,y)', hcat(x,y)',
                        [-1 -1], phs, pol, stc, K;
                        tree = tr, idx = ind_nn)[[1,2]]

# # These are only needed on the boundary nodes
# jj = jj[:,bb]
# wx = wx[:,bb]
# wy = wy[:,bb]
# awhv = awhv[:,bb]

###############################################################################

# Get initial conditions, which are also background states
pi_0, th_0, rho_0, e_0, q_0, p_0, T_0 = getInitialConditions(qq, Cv, R)

# Get initial averages for density, tracer density, and total energy
rho_0_avg = ((1/4)*sum(rho_0[cs1]) + (2/4)*sum(rho_0[cs2]) +
             (3/4)*sum(rho_0[cs3]) + (4/4)*sum(rho_0[cs4])) / numSquares
rhoq_0 = rho_0 .* q_0
rhoq_0_avg = ((1/4)*sum(rhoq_0[cs1]) + (2/4)*sum(rhoq_0[cs2]) +
              (3/4)*sum(rhoq_0[cs3]) + (4/4)*sum(rhoq_0[cs4])) / numSquares
E_0 = rho_0 .* e_0
E_0_avg = ((1/4)*sum(E_0[cs1]) + (2/4)*sum(E_0[cs2]) +
           (3/4)*sum(E_0[cs3]) + (4/4)*sum(E_0[cs4])) / numSquares

io = open("./results/pi_0.txt", "w")
writedlm(io, pi_0, ' ')
close(io)

io = open("./results/th_0.txt", "w")
writedlm(io, th_0, ' ')
close(io)

io = open("./results/rho_0.txt", "w")
writedlm(io, rho_0, ' ')
close(io)

io = open("./results/e_0.txt", "w")
writedlm(io, e_0, ' ')
close(io)

io = open("./results/p_0.txt", "w")
writedlm(io, p_0, ' ')
close(io)

io = open("./results/T_0.txt", "w")
writedlm(io, T_0, ' ')
close(io)

# Set initial conditions for main solution array U
U = zeros(length(x), 5)
U[:,1] = pi_0
U[:,2] .= 0
U[:,3] .= 0
U[:,4] = th_0
U[:,5] = q_0

###############################################################################

# Initialize dummy arrays to be used in Runge-Kutta or Adams-Bashforth

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
                                       Wx, Wy, aWhv, pi_0, th_0,
                                       Tx, Ty, bb, ff, fanSpeed,
                                       # TX, TY, wx, wy, awhv, jj,
                                       Cp, Cv, R)

###############################################################################

function printAndSave!(i, U, xt, yt, frame,
                       maxVel, minPi, maxPi, minU, maxU, minV, maxV,
                       minTh, maxTh, minQ, maxQ, minT, maxT)
    
    maxVel = vcat(maxVel, maximum(sqrt.(U[:,2] .^ 2 .+ U[:,3] .^ 2)))
    tmp = U[:,1] - pi_0
    minPi  = vcat(minPi,  minimum(tmp))
    maxPi  = vcat(maxPi,  maximum(tmp))
    minU   = vcat(minU,   minimum(U[:,2]))
    maxU   = vcat(maxU,   maximum(U[:,2]))
    minV   = vcat(minV,   minimum(U[:,3]))
    maxV   = vcat(maxV,   maximum(U[:,3]))
    tmp = U[:,4] - th_0
    minTh  = vcat(minTh,  minimum(tmp))
    maxTh  = vcat(maxTh,  maximum(tmp))
    minQ   = vcat(minQ,   minimum(U[:,5]))
    maxQ   = vcat(maxQ,   maximum(U[:,5]))
    tmp = U[:,1] .* U[:,4] .- T_0
    minT = vcat(minT, minimum(tmp))
    maxT = vcat(maxT, maximum(tmp))

    @printf("frame = %.0f,  t = %.5f\n", frame, t[i])
    @printf("maxPi = %.5f,  maxU = %.5f,  maxV = %.5f,  maxTh = %.5f,  maxQ = %.5f\n", 
            maxPi[end],
            maxU[end],
            maxV[end],
            maxTh[end],
            maxQ[end])
    @printf("minPi = %.5f,  minU = %.5f,  minV = %.5f,  minTh = %.5f,  minQ = %.5f\n\n", 
            minPi[end],
            minU[end],
            minV[end],
            minTh[end],
            minQ[end])

    io = open(@sprintf("./results/pi_%04d.txt",frame), "w")
    writedlm(io, U[:,1], ' ')
    close(io)

    io = open(@sprintf("./results/u_%04d.txt",frame), "w")
    writedlm(io, U[:,2], ' ')
    close(io)

    io = open(@sprintf("./results/v_%04d.txt",frame), "w")
    writedlm(io, U[:,3], ' ')
    close(io)
    
    io = open(@sprintf("./results/th_%04d.txt",frame), "w")
    writedlm(io, U[:,4], ' ')
    close(io)

    io = open(@sprintf("./results/q_%04d.txt",frame), "w")
    writedlm(io, U[:,5], ' ')
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
        maximum(abs.(U[:,2])) > 100 ||
        maximum(abs.(U[:,3])) > 100 ||
        maximum(abs.(U[:,4])) > 1000 ||
        !iszero(isnan.(U)))
        println("It blew up.")
        stop = true
    end

    return stop, maxVel, minPi, maxPi, minU, maxU, minV, maxV, minTh, maxTh,
           minQ, maxQ, minT, maxT

end

###############################################################################

# The main time-stepping loop

function mainLoop!(U, f1, f2, f3)

    # Some things to keep track of periodically
    maxVel = []
    minPi  = [];  maxPi = []
    minU   = [];  maxU  = []
    minV   = [];  maxV  = []
    minTh  = [];  maxTh = []
    minQ   = [];  maxQ  = []
    minT   = [];  maxT  = []

    # Save initial stuff
    frame = 0
    stop, maxVel, minPi, maxPi, minU, maxU, minV, maxV, minTh, maxTh, minQ,
    maxQ, minT, maxT = printAndSave!(1, U, x, y, frame, maxVel,
                                     minPi, maxPi, minU, maxU, minV, maxV,
                                     minTh, maxTh, minQ, maxQ, minT, maxT)
    frame = 1

    # Testing out some methods for keeping track of moving points.  (xt,yt)
    # are the locations of moving particles, which will follow the flow.
    DT = tf / (numSaves - 1)
    xt = copy(x)
    yt = copy(y)
    u0 = copy(U[:,2])
    v0 = copy(U[:,3])

    # Do two steps of Runge-Kutta for prog vars, to get first three time-levels
    # U1 = copy(U)
    # f1 = odefun!(t[1], U, f1)
    U = rk3!(t[1], U, odefun!, dt, q1, q2)
    # U2 = copy(U)
    # f2 = odefun!(t[2], U, f2)
    U = rk3!(t[2], U, odefun!, dt, q1, q2)
    # U3 = copy(U)
    # f3 = odefun!(t[3], U, f3)

    for i in 3 : Int(tf/dt) + 1

        if mod(i-1, Int(tf/dt/(numSaves-1))) == 0
            # Get the velocity at the half time level (used twice below)
            uHalf = (u0 .+ U[:,2]) ./ 2
            vHalf = (v0 .+ U[:,3]) ./ 2
            # Let the moving points take a leap forward to the present time
            W = getDM(hcat(xt,yt)', hcat(x,y)',
                      [0 0], 1, 0, stc, 0;
                      tree = tr)[1]
            xHalf = xt .+ (DT/2) .* (W * u0)
            yHalf = yt .+ (DT/2) .* (W * v0)
            W = getDM(hcat(xHalf,yHalf)', hcat(x,y)',
                      [0 0], 1, 0, stc, 0;
                      tree = tr)[1]
            xt = xt .+ DT .* (W * uHalf)
            yt = yt .+ DT .* (W * vHalf)
            # # Remove wayward points (this has been causing stationary pts)
            # ind_in = (xt .>= 1 - overlap) .& (xt .<= wid + overlap) .&
            #          (yt .>= 1 - overlap) .& (yt .<= len + overlap)
            # xt = copy(xt[ind_in])
            # yt = copy(yt[ind_in])
            # Do one step of semi-lagrangian tracer transport
            xHalf = x .+ (DT/2) .* u0
            yHalf = y .+ (DT/2) .* v0
            W = getDM(hcat(xHalf,yHalf)', hcat(x,y)',
                      [0 0], phs, pol, stc, 0;
                      tree = tr)[1]
            xq = x .+ DT .* (W * uHalf)
            yq = y .+ DT .* (W * vHalf)
            W = getDM(hcat(x,y)', hcat(xq,yq)',
                      [0 0], 1, 0, stc, 0)[1]
            U[:,5] = W * U[:,5]
            # Shift old velocity to the new time level
            u0 = U[:,2]
            v0 = U[:,3]
            # Apply a fixer for mass and tracer mass
            U = simpleFixer!(U, rho_0_avg, rhoq_0_avg, E_0_avg,
                             cs1, cs2, cs3, cs4, numSquares,
                             Cv, R, p_0)
            # Print some info and save some things
            stop, maxVel, minPi, maxPi, minU, maxU, minV, maxV, minTh, maxTh,
            minQ, maxQ, minT, maxT = printAndSave!(i, U, xt, yt, frame,
                                                   maxVel, minPi, maxPi, minU,
                                                   maxU, minV, maxV, minTh,
                                                   maxTh, minQ, maxQ, minT,
                                                   maxT)
            # Time the update
            @time begin
                U = rk3!(t[i], U, odefun!, dt, q1, q2)
            end
            # Stop if it blew up
            if stop
                break
            end
            frame = frame + 1
        else
            # Just step forward in time
            U = rk3!(t[i], U, odefun!, dt, q1, q2)
        end

    end
    
    io = open("./results/maxVel.txt", "w")
    writedlm(io, maximum(maxVel), ' ')
    close(io)

    io = open("./results/minPi.txt", "w")
    writedlm(io, minimum(minPi), ' ')
    close(io)

    io = open("./results/maxPi.txt", "w")
    writedlm(io, maximum(maxPi), ' ')
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

    io = open("./results/minTh.txt", "w")
    writedlm(io, minimum(minTh), ' ')
    close(io)

    io = open("./results/maxTh.txt", "w")
    writedlm(io, maximum(maxTh), ' ')
    close(io)

    io = open("./results/minQ.txt", "w")
    writedlm(io, minimum(minQ), ' ')
    close(io)

    io = open("./results/maxQ.txt", "w")
    writedlm(io, maximum(maxQ), ' ')
    close(io)

    io = open("./results/minT.txt", "w")
    writedlm(io, minimum(minT), ' ')
    close(io)

    io = open("./results/maxT.txt", "w")
    writedlm(io, maximum(maxT), ' ')
    close(io)

    return U

end

U = mainLoop!(U, f1, f2, f3)

###############################################################################













