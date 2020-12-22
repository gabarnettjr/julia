
# include("../disk/makeRadialNodes.jl")
# include("../disk/appendGhostNodes.jl")

function makeNodes!(layers, radius)

    x1, y1 = makeRadialNodes(layers - 2)
    dr = 1 / (layers - 3)
    x1, y1 = appendGhostNodes!(x1, y1, dr)
    x1 = radius .* x1
    y1 = radius .* y1
    dr = sqrt((x1[1] - x1[2])^2 + (y1[1] - y1[2])^2)
    radius = radius + 3 ./ 2 .* dr

    x2 = copy(x1) .- (2*radius + 3*dr)
    y2 = copy(y1)
    
    nb = (layers - 3) * 6
    y_connect = y1[Int(end - nb/2 - 2) : Int(end - nb/2 + 4)]
    # y_connect = (y_connect[1:end-1] .+ y_connect[2:end]) ./ 2
    x_connect = (-radius - dr/2 - dr) .* ones(size(y_connect))
    x1 = vcat(x1, x_connect)
    y1 = vcat(y1, y_connect)

    x3 = copy(x1) .+ (2*radius + 3*dr)
    y3 = copy(y1)

    x = vcat(x1, x2, x3)
    y = vcat(y1, y2, y3)

    return x, y, dr, radius

end

