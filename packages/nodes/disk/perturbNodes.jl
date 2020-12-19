
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
