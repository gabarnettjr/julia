
function getIndices(x, y, dr, s, doorWidth, outAngle, dth)

    bb = []
    ff = []
    ii = []

    for i in 1 : length(x)

        r = sqrt(x[i]^2 + y[i]^2)
        th = atan(y[i], x[i])

        if abs(r - radius) < dr/2

            bb = vcat(bb, i)

        elseif r > s*radius - dr/2 && r < s*radius + 4*dr + dr/2 &&
        abs(th) < pi - doorWidth/2 + dth/2 &&
        abs(th - outAngle) > doorWidth/2 - dth/2

            bb = vcat(bb, i)
        
        elseif abs(r - s*radius) < dr/2 &&
        abs(th) > pi - doorWidth/2 + dth/2
            
            ff = vcat(ff, i)

        else

            ii = vcat(ii, i)

        end

    end

    return bb, ff, ii

end

