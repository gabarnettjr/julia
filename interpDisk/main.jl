
using DelimitedFiles
using HaltonSequences

#####################################################################

include("../packages/phs2.jl")
include("../packages/disk/nodes.jl")

#####################################################################

frq = 3
n = 1000
layers = 65
phs = 5
pol = 2
stc = 19
del = .2

#####################################################################

# The Halton nodes where information is known
H = HaltonPoint(2, length=n)
x = zeros(n)
y = zeros(n)
for i in 1:n
    x[i] = H[i][1]
    y[i] = H[i][2]
end
x = -1-del/2 .+ (2+del) * x
y = -1-del/2 .+ (2+del) * y

# # The random nodes where information is known
# p = -1-del/2 .+ (2+del) * rand(n,2)
# x = p[:,1]
# y = p[:,2]

# # The Cartesian nodes where information is known
# x = range(-1-del/2, stop=1+del/2, length=Int(round(sqrt(n))))
# y = copy(x)
# xx = repeat(x, 1, length(y))
# yy = repeat(y', length(x), 1)
# x = xx[:]
# y = yy[:]

# Only keep nodes inside a certain radius
ii = sqrt.(x .^ 2 .+ y .^ 2) .<= 1+del/2
x = x[ii]
y = y[ii]

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














