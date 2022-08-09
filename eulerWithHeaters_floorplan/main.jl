
using Printf
using LinearAlgebra
using DelimitedFiles

###############################################################################

# USER INPUT

# Final time (power of 2)
const tf = 2^2

# The type of floorplan
const floorplan = "square"

# Approximate node spacing in meters is dr, which is 1/k.
# The wall width is always 1, but dr controls how many layers wide each wall is.
# dr=1/7 => 6 layers wide.
# dr should be one divided by an odd number (5, 7, 9, ...).
# This is so the middle of a given block is a single row or column.
const k = 11
const dr = 1/k
# const overlap = Int((k-1)/2) * dr

###############################################################################

include("../packages/nodes/getFloorplanMatrix.jl")
include("../packages/nodes/expandFloorplanMatrix.jl")
include("../packages/nodes/getFloorplanNodesAndIndices.jl")
include("../packages/nodes/findCornersAndSides.jl")
include("../packages/nodes/tangents.jl")
include("../packages/nodes/countIntegrationSquares.jl")

include("../packages/eulerEquations/getConstants.jl")
include("../packages/eulerEquations/fanPrime.jl")
include("../packages/eulerEquations/heatPrime.jl")

include("../packages/phs2.jl")

# A structure to hold the constant inputs passed to rk3 and ODEfunction
struct constantInputs
    dt::Float64                    # delta t
    p_0::Array{Float64,1}          # background pressure
    rho_0::Array{Float64,1}        # background density
    T_0::Array{Float64,1}          # background temperature
    Tx::Array{Float64,1}           # tangent vecs on bndry, x-comp
    Ty::Array{Float64,1}           # tangent vecs on bndry, y-comp
    bb::Array{UInt,1}              # index of bndry pts
    ff2::Array{UInt,1}             # fan index west-to-east
    ff3::Array{UInt,1}             # fan index south-to-north
    ff4::Array{UInt,1}             # fan index east-to-west
    ff5::Array{UInt,1}             # fan index north-to-south
    fanSpeed::Float64              # limiting wind speed at fan
    tempDiff::Float64              # diff between heater and background
    Cp::Float64                    # specific heat at constant pressure
    Cv::Float64                    # specific heat at constant volume
    R::Float64                     # gas constant for dry air
    jj::Array{UInt,2}              # indices of nearest neighbors
    wx::Array{Float64,2}           # weights for d/dx
    wy::Array{Float64,2}           # weights for d/dy
    awhv::Array{Float64,2}         # weights for hyperviscosity
end

###############################################################################

# Functions for manipulating objects of type q

# A mutable structure for holding the main solution and similar objects
mutable struct q
    u::Array{Float64,1}
    v::Array{Float64,1}
    rho::Array{Float64,1}
    T::Array{Float64,1}
end

# A function for adding objects of type q
function qAdd(U::q, V::q, W::q)
    W.u   .= U.u   .+ V.u
    W.v   .= U.v   .+ V.v
    W.rho .= U.rho .+ V.rho
    W.T   .= U.T   .+ V.T
    return W
end

# A function for multiplying objects of type q by a scalar
function qMult(r::Float64, U::q, W::q)
    W.u   .= r .* U.u
    W.v   .= r .* U.v
    W.rho .= r .* U.rho
    W.T   .= r .* U.T
    return W
end

###############################################################################

include("getInitialConditions.jl")
include("simpleFixer.jl")
include("ODEfunction.jl")
include("rk3.jl")

###############################################################################

# This block contains parameters that are seldom changed.

# Constant background temperature, in Kelvin
const bgTemp = 300.
# const bgTemp = 280.

# Temperature difference between heater (fan) and background state, in Kelvin:
const tempDiff = 0.
# const tempDiff = 20.

# Speed of air blown by the fan (m/s)
const fanSpeed = 10.

# Delta t in seconds (power of 2 times dr)
const dt = 2^-11 * dr

# Array of all times in seconds
const t = range(0, stop=tf, step=dt)

# The number of times to save (power of 2 PLUS ONE)
const numSaves = 2^(Int(log2(tf)) + 6) + 1

# Exponent in the polyharmonic spline function
const phs = 5

# Highest degree of polynomial to include in the basis
const pol = 2

# Stencil size (15 is largest you can get away with when k = 5)
const stc = 15

# Hyperviscosity exponent
const K = 2

# Set the hyperviscosity coefficient
if K == 2
    # const a = 0.
    const a = -2^3 * dr ^ (2*K-1)
elseif K == 3
    const a = 2^0 * dr ^ (2*K-1)
else
    error("K should be 2 or 3.")
end

###############################################################################

# Remove old files and remake the results directory
rm("../../../results", recursive = true)
mkdir("../../../results")

# Array of times when things are saved
const tSaves = range(0, stop = tf, length = numSaves)

# Larger delta-t for time-stepping the moving points (just for visualization)
const DT = tf / (numSaves - 1)

# Get the nodes on the domain, with the door built in
const M_0, len, wid = getFloorplanMatrix(plan = floorplan)
const M = expandFloorplanMatrix(M_0, Int(1/dr))
const x, y, xNon, yNon, bb, ff2, ff3, ff4, ff5, ii =
    getFloorplanNodesAndIndices(M, len, wid, dr)

# Find corner and side nodes to be used in mass fixer(s)
const cs1, cs2, cs3, cs4 = findCornersAndSides(x, y, dr)

# Count the number of small squares, to be used in the mass fixer(s)
const numSquares = countIntegrationSquares(x, y, len, wid, dr)

# Get the unit tangent vector at each boundary node
const Tx, Ty = tangents(x[bb], y[bb])

# Physical constants having to do with the fluid (dry air)
const Cp, Cv, R = getConstants()

###############################################################################

# Save things that will be needed in the python plotting routine.

function save(fileName, variable)
    io = open(fileName, "w")
    writedlm(io, variable, ' ')
    close(io)
end

save("../../../results/x.txt",           x)
save("../../../results/y.txt",           y)
save("../../../results/dr.txt",         dr)
save("../../../results/Cp.txt",         Cp)
save("../../../results/Cv.txt",         Cv)
save("../../../results/R.txt",           R)
save("../../../results/tSaves.txt", tSaves)
save("../../../results/bb.txt",         bb)
save("../../../results/ff2.txt",       ff2)
save("../../../results/ff3.txt",       ff3)
save("../../../results/ff4.txt",       ff4)
save("../../../results/ff5.txt",       ff5)
save("../../../results/xNon.txt",     xNon)
save("../../../results/yNon.txt",     yNon)

###############################################################################

# Get all of the weights that will be needed in the ODE function

const wx, jj, tr, ind_nn = phs2_getDM(hcat(x,y)', hcat(x,y)', [1 0], phs, pol,
                                      stc, 0)[2:5]

const wy = phs2_getDM(hcat(x,y)', hcat(x,y)', [0 1], phs, pol, stc, 0;
                      tree = tr, idx = ind_nn)[2]

const awhv = a .* phs2_getDM(hcat(x,y)', hcat(x,y)', [-1 -1], phs, pol, stc, K;
                             tree = tr, idx = ind_nn)[2]

###############################################################################

# Get initial conditions, which are also background states
const rho_0, p_0, T_0 = getInitialConditions(R, bgTemp)

# Get initial averages for density, tracer density, and total energy
const rho_0_avg = ((1/4)*sum(rho_0[cs1]) + (2/4)*sum(rho_0[cs2]) +
             (3/4)*sum(rho_0[cs3]) + (4/4)*sum(rho_0[cs4])) / numSquares
# rhoq_0 = rho_0 .* q_0
# rhoq_0_avg = ((1/4)*sum(rhoq_0[cs1]) + (2/4)*sum(rhoq_0[cs2]) +
#               (3/4)*sum(rhoq_0[cs3]) + (4/4)*sum(rhoq_0[cs4])) / numSquares
# E_0 = rho_0 .* e_0
# E_0_avg = ((1/4)*sum(E_0[cs1]) + (2/4)*sum(E_0[cs2]) +
#            (3/4)*sum(E_0[cs3]) + (4/4)*sum(E_0[cs4])) / numSquares

save("../../../results/rho_0.txt", rho_0)
save("../../../results/p_0.txt",     p_0)
save("../../../results/T_0.txt",     T_0)

# Set initial conditions for main solution U
U = q(zeros(length(x)), zeros(length(x)), copy(rho_0), copy(T_0))

# Temporary array to be used inside rk3:
Utmp = q(zeros(length(x)), zeros(length(x)), copy(rho_0), copy(T_0))

# Initialize arrays to store perturbations and divergence inside ODEfunction
pPrime   = zeros(length(x))
rhoPrime = zeros(length(x))
Tprime   = zeros(length(x))
D        = zeros(length(x))

# Temporary array to remove memory allocations for matrix multiplications
tmp1 = zeros(length(x))
tmp2 = zeros(length(x))

###############################################################################

# Initialize dummy arrays to be used for time-stepping

q1 = q(zeros(length(x)), zeros(length(x)), copy(rho_0), copy(T_0))
q2 = q(zeros(length(x)), zeros(length(x)), copy(rho_0), copy(T_0))
# q3 = zeros(size(U))                                           #needed for RK4
# q4 = zeros(size(U))                                           #needed for RK4

# f1 = zeros(size(U))                                           #needed for AB3
# f2 = zeros(size(U))                                           #needed for AB3
# f3 = zeros(size(U))                                           #needed for AB3

###############################################################################

# Create an instance of the struct constantInputs, to hold all these.

ci = constantInputs(dt, p_0, rho_0, T_0, Tx, Ty, bb, ff2, ff3, ff4, ff5, fanSpeed, tempDiff, Cp, Cv, R, jj, wx, wy, awhv)

###############################################################################

function printAndSave!(i, U, xt, yt, frame, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT)

    fp1 = frame + 1

    maxVel[fp1] = maximum(sqrt.(U.u .^ 2 .+ U.v .^ 2))
    tmp = U.rho * R .* U.T - p_0
    minP[fp1] = minimum(tmp)
    maxP[fp1] = maximum(tmp)
    tmp = U.rho - rho_0
    minRho[fp1] = minimum(tmp)
    maxRho[fp1] = maximum(tmp)
    minU[fp1] = minimum(U.u)
    maxU[fp1] = maximum(U.u)
    minV[fp1] = minimum(U.v)
    maxV[fp1] = maximum(U.v)
    tmp = U.T - T_0
    minT[fp1] = minimum(tmp)
    maxT[fp1] = maximum(tmp)

    @printf("frame = %.0f,  t = %.5f\n", frame, t[i])

    @printf("maxP = %6.0f,  maxRho = %8.5f,  maxU = %5.1f,  maxV = %5.1f,  maxT = %7.4f\n", maxP[fp1], maxRho[fp1], maxU[fp1], maxV[fp1], maxT[fp1])

    @printf("minP = %6.0f,  minRho = %8.5f,  minU = %5.1f,  minV = %5.1f,  minT = %7.4f\n\n", minP[fp1], minRho[fp1], minU[fp1], minV[fp1], minT[fp1])

    save(@sprintf("../../../results/u_%04d.txt",frame),     U.u)
    save(@sprintf("../../../results/v_%04d.txt",frame),     U.v)
    save(@sprintf("../../../results/rho_%04d.txt",frame), U.rho)
    save(@sprintf("../../../results/T_%04d.txt",frame),     U.T)
    save(@sprintf("../../../results/xt_%04d.txt",frame),     xt)
    save(@sprintf("../../../results/yt_%04d.txt",frame),     yt)

    stop = false

    # Stop running if the numerical solution blows up
    if (maximum(abs.(U.u)) > 500 || maximum(abs.(U.v)) > 500 ||
    maximum(abs.(U.rho)) > 10 || maximum(abs.(U.T)) > 500)
        println("It blew up.")
        stop = true
    end

    return stop, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT

end

###############################################################################

# The main time-stepping loop

function mainLoop(U::q)

    # Some things to keep track of periodically
    maxVel = zeros(numSaves);
    minP   = zeros(numSaves);  maxP   = zeros(numSaves);
    minRho = zeros(numSaves);  maxRho = zeros(numSaves);
    minU   = zeros(numSaves);  maxU   = zeros(numSaves);
    minV   = zeros(numSaves);  maxV   = zeros(numSaves);
    minT   = zeros(numSaves);  maxT   = zeros(numSaves);

    # Save initial stuff
    frame = 0
    stop, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT = printAndSave!(1, U, x, y, frame, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT)
    frame = 1

    # Testing out some methods for keeping track of moving points.  (xt,yt)
    # are the locations of moving particles, which will follow the flow.
    xt = copy(x)
    yt = copy(y)
    u0 = copy(U.u)
    v0 = copy(U.v)

    # Do two steps of Runge-Kutta for prog vars, to get first three time-levels
    # U1 = copy(U)
    # f1 = odefun!(t[1], U, f1)
    U = rk3!(t[1], U,
             q1, q2, Utmp, pPrime, rhoPrime, Tprime, D,
             tmp1, tmp2, ci)
    # U2 = copy(U)
    # f2 = odefun!(t[2], U, f2)
    U = rk3!(t[2], U,
             q1, q2, Utmp, pPrime, rhoPrime, Tprime, D,
             tmp1, tmp2, ci)
    # U3 = copy(U)
    # f3 = odefun!(t[3], U, f3)

    for i in 3 : Int(tf/dt) + 1

        if mod(i-1, Int(tf/dt/(numSaves-1))) == 0
            # Get the velocity at the half time level
            uHalf = (u0 .+ U.u) ./ 2
            vHalf = (v0 .+ U.v) ./ 2
            # Let the moving points take a leap forward to the present time
            W = phs2_getDM(hcat(xt,yt)', hcat(x,y)', [0 0], 1, 0, stc, 0; tree = tr)[1]
            xHalf = xt .+ (DT/2) .* (W * u0)
            yHalf = yt .+ (DT/2) .* (W * v0)
            W = phs2_getDM(hcat(xHalf,yHalf)', hcat(x,y)', [0 0], 1, 0, stc, 0; tree = tr)[1]
            xt .= xt .+ DT .* (W * uHalf)
            yt .= yt .+ DT .* (W * vHalf)
            # # Remove wayward points (this has been causing stationary pts)
            # ind_in = (xt .>= 1 - overlap) .& (xt .<= wid + overlap) .&
            #          (yt .>= 1 - overlap) .& (yt .<= len + overlap)
            # xt = copy(xt[ind_in])
            # yt = copy(yt[ind_in])
            # # Do one step of semi-lagrangian tracer transport
            # xHalf = x .+ (DT/2) .* u0
            # yHalf = y .+ (DT/2) .* v0
            # W = phs2_getDM(hcat(xHalf,yHalf)', hcat(x,y)', [0 0], phs, pol, stc, 0; tree = tr)[1]
            # xq = x .+ DT .* (W * uHalf)
            # yq = y .+ DT .* (W * vHalf)
            # W = phs2_getDM(hcat(x,y)', hcat(xq,yq)', [0 0], 1, 0, stc, 0)[1]
            # U[:,5] = W * U[:,5]
            # Shift old velocity to the new time level
            u0 = copy(U.u)
            v0 = copy(U.v)
            # Apply a fixer for mass
            U = simpleFixer!(U, rho_0_avg, cs1, cs2, cs3, cs4, numSquares, Cv, R, p_0)
            # Print some info and save some things
            stop, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT = printAndSave!(i, U, xt, yt, frame, maxVel, minP, maxP, minRho, maxRho, minU, maxU, minV, maxV, minT, maxT)
            # Stop if it blew up
            if stop
                break
            end
            # Time the update
            @time begin
                U = rk3!(t[i], U,
                         q1, q2, Utmp, pPrime, rhoPrime, Tprime, D,
                         tmp1, tmp2, ci)
            end
            frame = frame + 1
        else
            # Just step forward in time
            U = rk3!(t[i], U,
                     q1, q2, Utmp, pPrime, rhoPrime, Tprime, D,
                     tmp1, tmp2, ci)
        end

    end
    
    save("../../../results/maxVel.txt", maximum(maxVel))
    save("../../../results/minP.txt",     minimum(minP))
    save("../../../results/maxP.txt",     maximum(maxP))
    save("../../../results/minRho.txt", minimum(minRho))
    save("../../../results/maxRho.txt", maximum(maxRho))
    save("../../../results/minU.txt",     minimum(minU))
    save("../../../results/maxU.txt",     maximum(maxU))
    save("../../../results/minV.txt",     minimum(minV))
    save("../../../results/maxV.txt",     maximum(maxV))
    save("../../../results/minT.txt",     minimum(minT))
    save("../../../results/maxT.txt",     maximum(maxT))

end

mainLoop(U)

###############################################################################

