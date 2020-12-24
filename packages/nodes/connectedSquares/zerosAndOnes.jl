
function zerosAndOnes(;domain = 1)
    
    if domain == 1
        
        wid = 20
        len = 20

        z1 = zeros(len + 2, wid + 2)
        z1[2:end-1, 2:end-1] .= 1

        z1[10:13, 10:13] .= 0

    elseif domain == 2

        wid = 16
        len = 8

        z1 = zeros(len + 2, wid + 2)
        
        # left room
        z1[2:end-1, [2,3,4,5]] .= 1
        
        # left hallway
        z1[[3,4],[6,7]] .= 1

        # middle room
        z1[2:end-1, [8,9,10,11]] .= 1

        # right hallway
        z1[[7,8], [12,13]] .= 1

        # right room
        z1[2:end-1, [14,15,16,17]] .= 1

    end

    return z1, wid, len

end

