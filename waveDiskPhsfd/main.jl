
using Printf
using LinearAlgebra
using DelimitedFiles

#####################################################################

include("../packages/disk/nodes.jl")

include("../packages/disk/waveEquation.jl")

include("../packages/phs2.jl")

#####################################################################

# USER INPUT

# Switch to decide whether to plot the eigenvalues of the matrix:
eigenvalues = true

# Switch to use alternate ODE function
useAlternateODEfunction = true

# Wave speed
c = 1/8

# Number of layers of radial nodes on the unit disk (odd number)
layers = 17

# Set how much the interior nodes will be perturbed
ptb = .35

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
tf = 50

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

# Perturb the nodes by percentage ptb of node spacing
x, y = perturbNodes!(x, y, layers, ptb, bb)

# Number of nodes total
n = length(x)

# Number of interior nodes
ni = length(x[ii])

# Number of boundary nodes
nb = length(x[bb])

# Get all of the DMs that will be needed for the wave equation
cWx = c * getDM(hcat(x,y)', hcat(x,y)', [1 0], phs, pol, stc, 0)
cWy = c * getDM(hcat(x,y)', hcat(x,y)', [0 1], phs, pol, stc, 0)
aWhv = a * getDM(hcat(x,y)', hcat(x,y)', [-1 -1], phs, pol, stc, K)

#####################################################################

# Construct the matrix that would apply the entire ODE function
# in a single matrix-vector multiply, and then get its eigenvalues

null = sparse(zeros(length(x), length(x)))

A = hcat(   aWhv[ii,ii],    cWx[ii,:], cWy[ii,:])
A = vcat(A, hcat(cWx[:,ii], null,      null))
A = vcat(A, hcat(cWy[:,ii], null,      null))

if eigenvalues

    eigenvalues = eigvals(Matrix(A))

    io = open("./results/e_real.txt", "w")
    writedlm(io, real(eigenvalues), ' ')
    close(io)

    io = open("./results/e_imag.txt", "w")
    writedlm(io, imag(eigenvalues), ' ')
    close(io)

end

#####################################################################

# Initialize main solution array U
x0 = 0.1
y0 = 0.2
if useAlternateODEfunction
    U = zeros(ni+n+n)
    U[1:ni] = exp.(-20*((x[ii] .- x0) .^ 2 .+ (y[ii] .- y0) .^ 2))
else
    U = zeros(length(x), 3)
    U[:,1] = exp.(-20*((x .- x0) .^ 2 .+ (y .- y0) .^ 2))
end

# Initialize dummy arrays to be used in Runge-Kutta
q1 = zeros(size(U))                              #needed in all cases
q2 = zeros(size(U))                           #needed for rk3 and rk4
q3 = []
q4 = []
if rkstages == 4
    q3 = zeros(size(U))                               #needed for rk4
    q4 = zeros(size(U))                               #needed for rk4
end

#####################################################################

# Define an ODE function which can be passed to rk

if useAlternateODEfunction

    function odefun!(t, U, A, dUdt)

        dUdt = ODEfunction2!(t, U, A, dUdt)

        return dUdt

    end

else

    function odefun!(t, U, A, dUdt)
    
        # Things that are needed from outside the scope of the function
        global cWx, cWy, aWhv, bb, ii
    
        dUdt = ODEfunction!(t, U, dUdt, cWx, cWy, aWhv, bb, ii)
    
        return dUdt
    
    end

end

#####################################################################

# Define the Runge-Kutta function based on the number of stages

if rkstages == 3

    function rk!(t, U, odefun!, dt, A, q1, q2, q3, q4)
        return rk3!(t, U, odefun!, dt, A, q1, q2, q3, q4)
    end

elseif rkstages == 4

    function rk!(t, U, odefun!, dt, q1, q2, q3, q4)
        return rk4!(t, U, odefun!, dt, q1, q2, q3, q4)
    end

else

    error("rkstages should be 3 or 4 please.")

end

#####################################################################

function printAndSave(i, t_i, U, ni, nb, n, frame)

    if useAlternateODEfunction
    
        @printf("i = %.0f,  t = %.5f\n", i, t[i])
        @printf("maxRho = %.5f,  maxU = %.5f,  maxV = %.5f\n", 
                maximum(U[1:ni]), maximum(U[ni+1:ni+n]),
                maximum(U[ni+n+1:ni+2*n]))
        @printf("minRho = %.5f,  minU = %.5f,  minV = %.5f\n\n", 
                minimum(U[1:ni]), minimum(U[ni+1:ni+n]),
                minimum(U[ni+n+1:ni+2*n]))

        io = open(@sprintf("./results/rho_%04d.txt",frame), "w")
        writedlm(io, vcat(U[1:ni],zeros(nb)), ' ')
        close(io)

        io = open(@sprintf("./results/u_%04d.txt",frame), "w")
        writedlm(io, U[ni+1:ni+n], ' ')
        close(io)

        io = open(@sprintf("./results/v_%04d.txt",frame), "w")
        writedlm(io, U[ni+n+1:ni+2*n], ' ')
        close(io)

    else

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

    end

end

#####################################################################

# The main time-stepping loop

frame = 0

function mainTimeSteppingLoop(rk!, odefun!, rkstages, U, t, dt
    , q1, q2, q3, q4, frame)

    for i in 1 : Int(tf/dt) + 1
    
        # Every once in a while, print some info and save some things
    
        if mod(i-1, Int((layers-1)/2)) == 0

            printAndSave(i, t[i], U, ni, nb, n, frame)
    
            frame = frame + 1
    
        end
    
        # Update the array U to the next time level
        # if mod(i-1, Int((layers-1)/2)) == 0
        #     @time begin
        #         U = rk!(t[i], U, odefun!, dt, A, q1, q2, q3, q4)
        #     end
        # else
            U = rk!(t[i], U, odefun!, dt, A, q1, q2, q3, q4)
        # end
        
        # Stop running if the numerical solution blows up
        if maximum(abs.(U)) > 10
            println("It blew up.")
            break
        end
    
    end

end

mainTimeSteppingLoop(rk!, odefun!, rkstages, U, t, dt,
                     q1, q2, q3, q4, frame)

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











