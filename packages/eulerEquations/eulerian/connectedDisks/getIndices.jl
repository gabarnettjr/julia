
function getIndices(x, y, dr, radius, doorWidth)
    
    indA = []
    indB = []
    indC = []

    v1 = []
    v2 = []
    
    for i in 1 : length(x)

        # The middle room
        th = atan(y[i], x[i])
        r = sqrt(x[i] ^2 + y[i] ^2)
        if ((th >= doorWidth/2) && (th <= pi - doorWidth/2)) ||
        ((th <= -doorWidth/2) && (th >= -pi + doorWidth/2))
            if abs(r - (radius + dr/2)) < 1e-3
                indC = vcat(indC, i)
            elseif abs(r - (radius - dr/2)) < 1e-3
                indB = vcat(indB, i)
            elseif abs(r - (radius - 3/2*dr)) < 1e-3
                indA = vcat(indA, i)
            end
        end
        if abs(x[i] + (radius + dr + dr/2)) < 1e-3
            v1 = vcat(v1, i)
        end
        
        # The left room
        th = atan(y[i], x[i] + (2*radius+3*dr))
        r = sqrt((x[i] + (2*radius+3*dr)) ^2 + y[i] ^2)
        if (th >= doorWidth/2) || (th <= -doorWidth/2)
            if abs(r - (radius + dr/2)) < 1e-3
                indC = vcat(indC, i)
            elseif abs(r - (radius - dr/2)) < 1e-3
                indB = vcat(indB, i)
            elseif abs(r - (radius - 3/2*dr)) < 1e-3
                indA = vcat(indA, i)
            end
        end
        
        # The right room
        th = atan(y[i], x[i] - (2*radius+3*dr))
        r = sqrt((x[i] - (2*radius+3*dr)) ^2 + y[i] ^2)
        if (th >= -pi + doorWidth/2) && (th <= pi - doorWidth/2)
            if abs(r - (radius + dr/2)) < 1e-3
                indC = vcat(indC, i)
            elseif abs(r - (radius - dr/2)) < 1e-3
                indB = vcat(indB, i)
            elseif abs(r - (radius - 3/2*dr)) < 1e-3
                indA = vcat(indA, i)
            end
        end
        if abs(x[i] - (radius + dr + dr/2)) < 1e-3
            v2 = vcat(v2, i)
        end
        
    end

    indC = vcat(indC, v1[1], v1[7])
    indC = vcat(indC, v2[1], v2[7])

    indB = vcat(indB, v1[2], v1[6])
    indB = vcat(indB, v2[2], v2[6])

    indA = vcat(indA, v1[3], v1[5])
    indA = vcat(indA, v2[3], v2[5])
        
    return indA, indB, indC

end

