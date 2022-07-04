
# Gets the unit tangent vector at each boundary node (xb,yb).
# It works by finding the two nearest neighbors to each node and returning
# the vector that connects these two points, normalized to  magnitude 1.

function tangents(xb::Array{Float64,1}, yb::Array{Float64,1})

    tree = KDTree(hcat(xb,yb)')

    idx = knn(tree, hcat(xb,yb)', 3, true)[1]

    Tx = zeros(size(xb))
    Ty = zeros(size(xb))

    for i in 1 : length(xb)
        Tx[i] = xb[idx[i][3]] - xb[idx[i][2]]
        Ty[i] = yb[idx[i][3]] - yb[idx[i][2]]
    end

    mag = sqrt.(Tx .^2 .+ Ty .^ 2)

    return Tx./mag, Ty./mag

end

