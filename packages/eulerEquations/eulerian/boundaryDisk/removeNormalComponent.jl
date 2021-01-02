
function removeNormalComponent!(u, v, Tx, Ty)

    uDotT = u .* Tx .+ v .* Ty

    u = uDotT .* Tx
    v = uDotT .* Ty

    return u, v

end

