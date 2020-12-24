
function makeNodes(z1, refinement, wid, len)

    dx = 1 / 2 ^ (refinement - 1)

    m, n = size(z1)

    tmp = dx/2 + (2^(refinement-1) - 1) * dx

    x = range(-tmp, stop = wid + tmp, step = dx)
    y = range(-tmp, stop = len + tmp, step = dx)

    X = zeros(Float64, m, n)
    Y = zeros(Float64, m, n)

    for j in 1 : n
        for i in 1 : m
            X[i,j] = x[j]
            Y[i,j] = y[i]
        end
    end

    return X, Y, dx

end

