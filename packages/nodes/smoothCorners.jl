
# Given boundary nodes in rectangular geometry, this function finds the corner
# nodes and replaces them with a slightly different point, so that the corners
# will be "smoothed".  Each corner point is replaced by the average of its
# three nearest boundary neighbors (including itself).

function smoothCorners!(xb, yb)

    tree = KDTree(hcat(xb,yb)')

    idx = knn(tree, hcat(xb,yb)', 3, true)[1]

    xnew = zeros(size(xb))
    ynew = zeros(size(yb))

    for i in 1 : length(xb)
        Tx = xb[idx[i][3]] - xb[idx[i][2]]
        Ty = yb[idx[i][3]] - yb[idx[i][2]]
        if abs(Tx) > 1e-3 && abs(Ty) > 1e-3
            xnew[i] = (xb[idx[i][1]] + xb[idx[i][2]] + xb[idx[i][3]]) / 3
            ynew[i] = (yb[idx[i][1]] + yb[idx[i][2]] + yb[idx[i][3]]) / 3
        else
            xnew[i] = xb[idx[i][1]]
            ynew[i] = yb[idx[i][1]]
        end
    end

    return xnew, ynew

end

