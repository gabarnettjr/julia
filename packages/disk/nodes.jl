#####################################################################

using DelimitedFiles

#####################################################################

"""
Distribute nodes somewhat uniformly on the unit disk by starting in
the middle and working your way out in layers.
"""

function makeRadialNodes(layers; n=6)

    h = 1 / (layers - 1)

    p = [0 0]

    k = 0

    for i in 1 : layers-1
        k = k + n
        th = range(0, stop=2*pi, length=k+1)
        th = th[1:end-1]
        p = vcat(p, i * h * hcat(cos.(th), sin.(th)))
    end

    # io = open("layers.txt", "w")
    # writedlm(io, layers, ' ')
    # close(io)

    # io = open("nPts.txt", "w")
    # writedlm(io, size(p,1), ' ')
    # close(io)

    # io = open("x.txt", "w")
    # writedlm(io, p[:,1], ' ')
    # close(io)

    # io = open("y.txt", "w")
    # writedlm(io, p[:,2], ' ')
    # close(io)

    return p[:,1], p[:,2]

end

#####################################################################

function perturbInterior!(x, y, layers, ptb, ii)

    r = -1 .+ 2 * rand(length(x[ii]), 2)

    x[ii] = x[ii] .+ 1/(layers-1) * ptb * r[:,1]
    y[ii] = y[ii] .+ 1/(layers-1) * ptb * r[:,2]

    return x, y 
end

#####################################################################
