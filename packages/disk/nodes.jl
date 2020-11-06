# module nodes

# export makeRadialNodes, perturbNodes!

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

    return p[:,1], p[:,2]

end

#####################################################################

function perturbNodes!(x, y, layers, ptb, bb)

    ran = -1 .+ 2 * rand(length(x), 2)
    x = x .+ 1/(layers-1) * ptb * ran[:,1]
    y = y .+ 1/(layers-1) * ptb * ran[:,2]

    r = sqrt.( x[bb] .^ 2 + y[bb] .^ 2)
    x[bb] = x[bb] ./ r
    y[bb] = y[bb] ./ r

    # th = atan.(y[bb], x[bb])
    # ran2 = -1 .+ 2 * rand(size(th))
    # dth = 2 * pi / length(th)
    # x[bb] = copy(cos.(th .+ dth*ptb*ran2))
    # y[bb] = copy(sin.(th .+ dth*ptb*ran2))

    # Why the fuck is this not working?  God damn it, I just want
    # to change the angle of the boundary nodes by a random
    # percentage of delta theta.  It shouldn't be this hard.

    return x, y 

end

#####################################################################

# end
