
function removeNormalComponent!(dudt::Array{Float64,1}, dvdt::Array{Float64,1},
                                bb::Array{UInt,1}, Tx::Array{Float64,1},
                                Ty::Array{Float64,1}, uDotT::Array{Float64,1})

    uDotT .= dudt[bb] .* Tx .+ dvdt[bb] .* Ty

    dudt[bb] .= uDotT .* Tx
    dvdt[bb] .= uDotT .* Ty

    return dudt, dvdt

end

