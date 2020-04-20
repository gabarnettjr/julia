
using Printf
using LinearAlgebra

#####################################################################

include("../packages/disk/nodes.jl")

include("../packages/disk/waveEquation.jl")

include("../packages/rk.jl")

#####################################################################

# USER INPUT

# Switch to decide whether to plot the eigenvalues of the matrix:
eigenvalues = false

# Wave speed
c = 1/8

# Number of layers of radial nodes on the unit disk
layers = 65

# Delta t
dt = 1/4/c * 1/(layers - 1)

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
tf = 100

#####################################################################

# Remove old files and remake the results directory
rm("results", recursive = true)
mkdir("results")

# Array of all times
t = range(0, stop=tf, step=dt)

# Get the nodes on the unit disk
x, y = makeRadialNodes(layers)

# Set the hyperviscosity coefficient
if K == 2
    a = -2^(-7) * c * 1 / (layers - 1) ^ (2*K-1)
elseif K == 3
    a = 2^(-10) * c * 1 / (layers - 1) ^ (2*K-1)
else
    error("Still need to implement other cases.")
end

# Index of the boundary nodes
bb = abs.(x .^ 2 .+ y .^ 2) .> (1 - 1e-6)

# Index of the interior nodes
ii = abs.(x .^ 2 .+ y .^ 2) .< (1 - 1e-6)

# Get all of the DMs that will be needed for the wave equation
cWx, cWy, aWhv = getAllDMs(c, hcat(x,y)', phs, pol, stc, K, a)

#####################################################################

# Construct the matrix that would apply the entier ODE function
# in a single matrix-vector multiply, and then get its eigenvalues

if eigenvalues

    null = zeros(length(x), length(x))

    A = hcat(Matrix(aWhv[ii,ii]), Matrix(cWx[ii,:]), Matrix(cWy[ii,:]))
    A = vcat(A, hcat(Matrix(cWx[:,ii]), null, null))
    A = vcat(A, hcat(Matrix(cWy[:,ii]), null, null))

    eigenvalues = eigvals(A)

    io = open("./results/e_real.txt", "w")
    writedlm(io, real(eigenvalues), ' ')
    close(io)

    io = open("./results/e_imag.txt", "w")
    writedlm(io, imag(eigenvalues), ' ')
    close(io)

end

#####################################################################

# Initialize main solution array U
U = zeros(length(x), 3)
U = initialCondition!(x, y, U)

# Initialize dummy arrays to be used in Runge-Kutta
q1 = zeros(size(U))                              #needed in all cases
if rkstages >= 3
    q2 = zeros(size(U))                       #needed for rk3 and rk4
end
if rkstages == 4
    q3 = zeros(size(U))                               #needed for rk4
    q4 = zeros(size(U))                               #needed for rk4
end

#####################################################################

# Define an ODE function which can be passed to rk

function odefun!(t, U, dUdt)

    # Things that are needed from outside the scope of the function
    global cWx, cWy, aWhv, bb, ii

    dUdt = ODEfunction!(t, U, dUdt, cWx, cWy, aWhv, bb, ii)

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

frame = 0

for i in 1 : Int(tf/dt)+1

    # Things needed from outside the scope of the for loop
    global rk!, odefun!, rkstages, U, t, dt, q1, q2, q3, q4, frame

    # Every once in a while, print some info and save some things

    if mod(i-1, Int((layers-1)/2)) == 0

        @printf("i = %.0f,  t = %.5f\n", i, t[i])
        @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f\n", 
                maximum(U[:,1]), maximum(U[:,2]), maximum(U[:,3]))
        @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f\n\n", 
                minimum(U[:,1]), minimum(U[:,2]), minimum(U[:,3]))

        io = open(@sprintf("./results/rho_%04d.txt",frame), "w")
        writedlm(io, U[:,1], ' ')
        close(io)

        io = open(@sprintf("./results/u_%04d.txt",frame), "w")
        writedlm(io, U[:,2], ' ')
        close(io)

        io = open(@sprintf("./results/v_%04d.txt",frame), "w")
        writedlm(io, U[:,3], ' ')
        close(io)

        frame = frame + 1

    end

    # Update the array U to the next time level
    U = rk!(t[i], U, odefun!, dt)

    # Stop running if the numerical solution blows up
    if maximum(abs.(U)) > 20
        println("It blew up.")
        break
    end

end

#####################################################################

# Save things that will be needed in the python plotting routine.
# Later on there might be more added to this list, so that all of
# the important details will be known to python.

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

#####################################################################











