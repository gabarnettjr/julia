
# Supply the boundary nodes as input for the round room problem, and this
# should give back unit tangent and unit normal vectors at each of those
# boundary nodes.

function tangentsNormals(x, y, radius, dr)

    Tx = zeros(length(x))
    Ty = zeros(length(x))
    Nx = zeros(length(x))
    Ny = zeros(length(x))

    for i in 1 : length(x)
        if (abs(x[i] + (radius + 3/2*dr)) < 1e-3) ||
           (abs(x[i] - (radius + 3/2*dr)) < 1e-3)
            if y[i] > 0
                Tx[i] = -1
                Ty[i] = 0
                Nx[i] = 0
                Ny[i] = -1
            else
                Tx[i] = 1
                Ty[i] = 0
                Nx[i] = 0
                Ny[i] = 1
            end
        else
            th = atan(y[i], x[i])
            Nx[i] = cos(th)
            Ny[i] = sin(th)
            Tx[i] = Ny[i]
            Ty[i] = -Nx[i]
        end
    end

    return Tx, Ty, Nx, Ny

end

