
# Gets the unit tangent vector at each boundary node (x,y).
# It works by finding the two nearest neighbors to each node and returning
# the vector that connects these two points, normalized to  magnitude 1.

function tangents(x, y)

    tree = KDTree(hcat(x,y)')

    idx = knn(tree, hcat(x,y)', 3)[1]

    Tx = zeros(size(x))
    Ty = zeros(size(x))

    for i in 1 : length(x)
        Tx[i] = x[idx[i][3]] - x[idx[i][2]]
        Ty[i] = y[idx[i][3]] - y[idx[i][2]]
    end

    mag = sqrt.(Tx .^2 .+ Ty .^ 2)

    return Tx./mag, Ty./mag

end

