
function zerosAndOnes(;domain = 1)
    
    if domain == 1

        A = zeros(10, 18)
        
        # left room
        A[2:end-1, [2,3,4,5]] .= 1
        
        # left hallway
        A[[3,4],[6,7]] .= 1

        # middle room
        A[2:end-1, [8,9,10,11]] .= 1

        # right hallway
        A[[7,8], [12,13]] .= 1

        # right room
        A[2:end-1, [14,15,16,17]] .= 1

    end

    return A

end

