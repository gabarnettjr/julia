
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
