
# Since the nodes are not at the center of the square cells, we need to know
# how many square cells there are in the domain, so that we can divide by this
# number in the mass and energy fixers when computing average values.

function countIntegrationSquares(x::Array{Float64,1}, y::Array{Float64,1},
                                 len::Int64, wid::Int64, dr::Float64)

    k = Int64(1/dr)

    # Get the points in the middle of each small square
    xcv = range(1-Int((k-1)/2)*dr+dr/2, wid+Int((k-1)/2)*dr-dr/2, step = dr)
    ycv = range(1-Int((k-1)/2)*dr+dr/2, len+Int((k-1)/2)*dr-dr/2, step = dr)
    xxc = zeros(length(ycv), length(xcv))
    yyc = zeros(length(ycv), length(xcv))
    for j in 1 : length(xcv)
        for i in 1 : length(ycv)
            xxc[i,j] = xcv[j]
            yyc[i,j] = ycv[i]
        end
    end

    # Find the four nearest neighbors to each of these central points
    tree = KDTree(hcat(x,y)')
    dis = knn(tree, hcat(xxc[:],yyc[:])', 4, true)[2]

    # If the fourth neighbor is within dr of the central point, then count it,
    # because that means the square is in the interior of the domain.
    numSquares = 0
    for i in 1 : length(xxc[:])
        if dis[i][4] < dr
            numSquares = numSquares + 1
        end
    end

    return numSquares

end

