
function removeNormalComponent!(u::Array{Float64,1}, v::Array{Float64,1},
                                Tx::Array{Float64,1}, Ty::Array{Float64,1}, uDotT::Array{Float64,1})

    uDotT .= u .* Tx .+ v .* Ty

    u .= uDotT .* Tx
    v .= uDotT .* Ty

    return u, v

end

