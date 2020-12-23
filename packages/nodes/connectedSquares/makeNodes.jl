
function makeNodes(z1, refinement)

    dx = 1 / refinement

    m, n = size(z1)

    x = range(-dx/2, stop = -dx/2 + (n-1)*dx, step = dx)
    y = range(-dx/2, stop = -dx/2 + (m-1)*dx, step = dx)

    X = zeros(m, n)
    Y = zeros(m, n)

    for j in 1 : n
        for i in 1 : m
            X[i,j] = x[j]
            Y[i,j] = y[i]
        end
    end

    return X, Y, dx

end

