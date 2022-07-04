
# Based on the expanded floorplan matrix, this function figures out where the
# nodes are, and it forms the index-sets to be used later in the code to grab,
# for example, just the boundary nodes (bb).

function getFloorplanNodesAndIndices(Mexp::Array{Int,2}, len::Int,
                                     wid::Int, dr::Float64)

    (m, n) = size(Mexp)

    k = Int(1/dr)

    xv = range(1-Int((k-1)/2)*dr, stop = wid+Int((k-1)/2)*dr, length = n)
    yv = range(1-Int((k-1)/2)*dr, stop = len+Int((k-1)/2)*dr, length = m)
    xx = zeros(m, n)
    yy = zeros(m, n)
    for j in 1 : n
        for i in 1 : m
            xx[i,j] = xv[j]
            yy[i,j] = yv[i]
        end
    end

    x::Array{Float64,1} = []
    y::Array{Float64,1} = []
    xNon::Array{Float64,1} = []
    yNon::Array{Float64,1} = []
    bb::Array{Int,1} = []
    ii::Array{Int,1} = []
    ff2::Array{Int,1} = []
    ff3::Array{Int,1} = []
    ff4::Array{Int,1} = []
    ff5::Array{Int,1} = []

    k = 1

    for j in 1 : n
        for i in 1 : m
            if Mexp[i,j] == -1
                xNon = vcat(xNon, xx[i,j])
                yNon = vcat(yNon, yy[i,j])
            else
                x = vcat(x, xx[i,j])
                y = vcat(y, yy[i,j])
                if Mexp[i,j] == 0
                    # Boundary index
                    bb = vcat(bb, k)
                elseif Mexp[i,j] == 1
                    # Interior index
                    ii = vcat(ii, k)
                elseif Mexp[i,j] == 2
                    # Fan blowing east
                    ff2 = vcat(ff2, k)
                elseif Mexp[i,j] == 3
                    # Fan blowing north
                    ff3 = vcat(ff3, k)
                elseif Mexp[i,j] == 4
                    # Fan blowing west
                    ff4 = vcat(ff4, k)
                elseif Mexp[i,j] == 5
                    # Fan blowing south
                    ff5 = vcat(ff5, k)
                else
                    error("Some unusable number came up in the matrix.")
                end
                k = k + 1
            end
        end
    end

    # # Find nearest interior node to each boundary node, which will be used in
    # # the enforcement of a Neumann boundary condition for temperature.
    # tree = KDTree(transpose(hcat(x[ii], y[ii])))
    # bc = knn(tree, transpose(hcat(x[bb], y[bb])), 1, true)[1]
    # tmp = zeros(Int, length(bb), 1)
    # for i in 1 : length(bb)
    #     tmp[i] = bc[i][1]
    # end
    # bc = copy(tmp)

    return x, y, xNon, yNon, bb, ff2, ff3, ff4, ff5, ii

end
