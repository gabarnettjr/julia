
using DelimitedFiles
using HaltonSequences

#####################################################################

include("../packages/phs2.jl")
# using phs2
include("../packages/disk/nodes.jl")

#####################################################################

# Frequency of the exact solution
frq = 3

# Number of nodes
n = 250*4^1

# number of layers in the regular evaluation points
layers = 129

phs = 5
pol = 2
stc = 19

# How much to go beyond a radius of 1
del = .1

#####################################################################

# x, y, iterations = makeHaltonNodes(n)
x, y, iterations = makeRandomNodes(n)
# x, y, iterations = makeCartesianNodes(n)

println(iterations)

# Scale the nodes out over the boundary
x = (1 + del) * x
y = (1 + del) * y

# The regular nodes where information is desired
xe, ye = makeRadialNodes(layers)

# The weights matrix for interpolating to the regular nodes
@time begin
    W = getDM(hcat(xe,ye)', hcat(x,y)', [0 0], phs, pol, stc, 0)
end

# Test function
f(x,y) = cos.(frq*pi*x) .* sin.(frq*pi*y)

# Approximation on evaluation nodes
app = W * f(x,y)

# Exact answer on evaluation nodes
exact = f(xe,ye)

#####################################################################

io = open("./results/x.txt", "w")
writedlm(io, x, ' ')
close(io)

io = open("./results/y.txt", "w")
writedlm(io, y, ' ')
close(io)

io = open("./results/xe.txt", "w")
writedlm(io, xe, ' ')
close(io)

io = open("./results/ye.txt", "w")
writedlm(io, ye, ' ')
close(io)

io = open("./results/app.txt", "w")
writedlm(io, app, ' ')
close(io)

io = open("./results/exact.txt", "w")
writedlm(io, exact, ' ')
close(io)

#####################################################################

