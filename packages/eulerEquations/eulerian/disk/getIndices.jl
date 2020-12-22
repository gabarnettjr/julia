
function getIndices(x, y, dr, radius)
    
    indA = []
    indB = []
    indC = []
    
    for i in 1 : length(x)
    
        r = sqrt(x[i] ^2 + y[i] ^2)
        
        if abs(r - (radius + dr/2)) < 1e-3
            indC = vcat(indC, i)
        elseif abs(r - (radius - dr/2)) < 1e-3
            indB = vcat(indB, i)
        elseif abs(r - (radius - 3/2*dr)) < 1e-3
            indA = vcat(indA, i)
        end
        
    end
        
    return indA, indB, indC

end

