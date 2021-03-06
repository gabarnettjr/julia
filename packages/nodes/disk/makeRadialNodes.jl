"""
Distribute nodes somewhat uniformly on the unit disk by starting in
the middle and working your way out in layers.
"""
function makeRadialNodes(layers; n=6)

    h = 1 / (layers - 1)

    p = [0 0]

    k = 0

    for i in 1 : layers-1
        k = k + n
        th = range(0, stop=2*pi, length=k+1)
        th = th[1:end-1]
        p = vcat(p, i * h * hcat(cos.(th), sin.(th)))
    end

    return p[:,1], p[:,2]

end

