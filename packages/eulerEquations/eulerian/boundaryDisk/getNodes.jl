
function getNodes(layers, radius;
                  n=6,
                  s = 1/2,
                  doorWidth = (layers-1)/4,
                  outAngle = (layers-1)*3/4)

    dth = 2 * pi / (s*(layers - 1) * 6)
    doorWidth = doorWidth * dth
    outAngle = outAngle * dth

    h = 1 / (layers - 1)

    p = [0 0]
    p_non = [NaN NaN]

    k = 0

    for i in 1 : layers-1
        global tmp
        k = k + n
        if i <= s * (layers - 1) || i >= s * (layers - 1) + 5
            th = range(-pi, stop=pi, length=k+1)
            th = th[1:end-1]
            if i == (s * (layers - 1))
                tmp = copy(th)
            end
            th_non = []
        elseif i >= s * (layers - 1) + 1 && i <= s * (layers - 1) + 3
            th1_non = tmp[(tmp .> -pi + doorWidth/2 + dth/2) .&
                          (tmp .< outAngle - doorWidth/2 - dth/2)]
            th2_non = tmp[(tmp .> outAngle + doorWidth/2 + dth/2) .&
                          (tmp .< pi - doorWidth/2 - dth + dth/2)]
            th_non = vcat(th1_non, th2_non)
            th1 = tmp[(tmp .> -pi - dth/2) .&
                      (tmp .< -pi + doorWidth/2 + dth/2)]
            th2 = tmp[(tmp .> outAngle - doorWidth/2 - dth/2) .&
                      (tmp .< outAngle + doorWidth/2 + dth/2)]
            th3 = tmp[(tmp .> pi - doorWidth/2 - dth/2) .&
                      (tmp .< pi - dth/2)]
            th = vcat(th1, th2, th3)
        elseif i == s * (layers - 1) + 4
            th = copy(tmp)
            th_non = []
        end
        p = vcat(p, i * h * hcat(cos.(th), sin.(th)))
        p_non = vcat(p_non, i * h * hcat(cos.(th_non), sin.(th_non)))
    end

    p = p .* radius
    p_non = p_non .* radius

    return p[:,1], p[:,2], p_non[2:end,1], p_non[2:end,2],
           s, doorWidth, outAngle, dth

end

