
function appendGhostNodes!(x, y, dr)

    xb = []
    yb = []

    for i in 1 : length(x)
        if abs(sqrt(x[i]^2 + y[i]^2) - 1) < 1e-3
            xb = vcat(xb, x[i])
            yb = vcat(yb, y[i])
        end
    end

    th = atan.(yb, xb)

    x = vcat(x, xb .+ dr .* cos.(th))
    y = vcat(y, yb .+ dr .* sin.(th))

    x = vcat(x, xb .+ 2 .* dr .* cos.(th))
    y = vcat(y, yb .+ 2 .* dr .* sin.(th))

    return x, y

end

