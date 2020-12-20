
# You give the coordinates of the smooth boundary points, in counterclockwise
# order, and this will give you back approximations of the unit tangent and
# unit normal vectors at each of those boundary points.

function tangentsNormals(x, y)

    Tx = zeros(length(x))
    Ty = zeros(length(x))
    Nx = zeros(length(x))
    Ny = zeros(length(x))

    Tx[1] = x[2] - x[end]
    Ty[1] = y[2] - y[end]

    for i in 2 : length(x) - 1
        Tx[i] = x[i+1] - x[i-1]
        Ty[i] = y[i+1] - y[i-1]
    end

    Tx[end] = x[1] - x[end-1]
    Ty[end] = y[1] - y[end-1]

    r = sqrt.(Tx .^2 .+ Ty .^2)
    Tx = Tx ./ r
    Ty = Ty ./ r

    Nx = -Ty
    Ny = Tx

    return Tx, Ty, Nx, Ny

end

