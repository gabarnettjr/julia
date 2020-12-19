
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
