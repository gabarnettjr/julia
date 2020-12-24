
function vecsToMats!(U, bool_all, nanMats)

    for i in 1 : 4
        tmp = nanMats[:,:,i]
        tmp[bool_all] = U[:,i]
        nanMats[:,:,i] = tmp
    end

    return nanMats

end

