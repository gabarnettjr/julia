
function matsToVecs!(U, bool_all, zeroVecs)

    for i in 1 : 4
        tmp = U[:,:,i][bool_all]
        zeroVecs[:,i] = tmp
    end

    return zeroVecs

end

