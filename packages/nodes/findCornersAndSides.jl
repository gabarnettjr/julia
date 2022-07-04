
# This function determines how many small squares each node is a part of.
# The number of squares that a given node is a part of determines how heavily
# weighted that node will be when estimating the average value of a function.
# This information is only ever used in simpleFixer.jl.

function findCornersAndSides(x, y, dr)

    tree = KDTree(hcat(x,y)')

    dis = knn(tree, hcat(x,y)', 9, true)[2]

    # cs stands for count squares
    cs1 = []
    cs2 = []
    cs3 = []
    cs4 = []

    for i in 1 : length(x)
        if (dis[i][4] < (3/2)*dr) && (dis[i][5] > (3/2)*dr)
            # node is part of only one square (corner node)
            cs1 = vcat(cs1, i)
        elseif (dis[i][6] < (3/2)*dr) && (dis[i][7] > (3/2)*dr)
            # node is part of two squares (side node)
            cs2 = vcat(cs2, i)
        elseif (dis[i][8] < (3/2)*dr) && (dis[i][9] > (3/2)*dr)
            # node is part of three squares (corner node)
            cs3 = vcat(cs3, i)
        else
            # most nodes are part of four squares (interior node)
            cs4 = vcat(cs4, i)
        end
    end

    return cs1, cs2, cs3, cs4

end
