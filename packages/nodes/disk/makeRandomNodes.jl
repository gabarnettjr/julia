
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
