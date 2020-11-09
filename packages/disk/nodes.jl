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

function makeCartesianNodes(n)

    x = []
    y = []
    j = Int(round(sqrt(n)))

    while length(x) < n
        x = range(-1, stop=1, length=j)
        y = copy(x)
        xx = repeat(x, 1, length(y))
        yy = repeat(y', length(x), 1)
        x = xx[:]
        y = yy[:]
        ii = sqrt.(x .^ 2 .+ y .^ 2) .<= 1
        x = x[ii]
        y = y[ii]
        j = j + 1
    end

    return x, y, j - Int(round(sqrt(n)))

end

#####################################################################

function makeRandomNodes(n)

    x = []
    y = []
    iterations = 0

    while length(x) != n
        p = -1 .+ 2 * rand(Int(round(4/pi*n)), 2)
        x = p[:,1]
        y = p[:,2]
        ii = sqrt.(x .^ 2 .+ y .^ 2) .<= 1
        x = x[ii]
        y = y[ii]
        iterations = iterations + 1
    end

    return x, y, iterations

end

#####################################################################

function makeHaltonNodes(n)

    x = []
    y = []
    ell = n

    while length(x) != n
        H = HaltonPoint(2, length=ell)
        x = zeros(ell)
        y = zeros(ell)
        for i in 1:ell
            x[i] = H[i][1]
            y[i] = H[i][2]
        end
        x = -1 .+ 2 * x
        y = -1 .+ 2 * y
        ii = sqrt.(x .^ 2 .+ y .^ 2) .<= 1
        x = x[ii]
        y = y[ii]
        ell = ell + 1;
    end

    return x, y, ell - n

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
