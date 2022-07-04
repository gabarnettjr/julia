
# include("removeNormalComponent.jl")
# include("fanPrime.jl")

function ODEfunction!(t, U, dUdt,
                      Wx, Wy, aWhv, pi_0, th_0,
                      Tx, Ty, bb, ff, fanSpeed,
                      # TX, TY, wx, wy, awhv, jj,
                      Cp, Cv, R)

    # (d/dt)(pi)
    pi_0 = U[:,1] .- pi_0
    dUdt[:,1] = -U[:,2] .* (Wx * pi_0) .-
                 U[:,3] .* (Wy * pi_0) .-
                 (R/Cv) .* U[:,1] .* (Wx * U[:,2] .+ Wy * U[:,3]) .+
                 aWhv * pi_0

    # (d/dt)(u)
    dUdt[:,2] = -U[:,2] .* (Wx * U[:,2]) .- U[:,3] .* (Wy * U[:,2]) .-
                 Cp .* U[:,4] .* (Wx * pi_0) .+
                 aWhv * U[:,2]

    # (d/dt)(v)
    dUdt[:,3] = -U[:,2] .* (Wx * U[:,3]) .- U[:,3] .* (Wy * U[:,3]) .-
                 Cp .* U[:,4] .* (Wy * pi_0) .+
                 aWhv * U[:,3]
    
    # # (d/dt)(th)
    # th_0 = U[:,4] .- th_0
    # dUdt[:,4] = -U[:,2] .* (Wx * th_0) .-
    #              U[:,3] .* (Wy * th_0) .+
    #              aWhv * th_0

    # # (d/dt)(q)
    # dUdt[:,5] = -U[:,2] .* (Wx * U[:,5]) .-
    #              U[:,3] .* (Wy * U[:,5]) .+
    #              aWhv * U[:,5]

    # Internal forcing (fans blowing air around the rooms)
    dUdt[ff,2] .= fanPrime(t; fs = fanSpeed)
    dUdt[ff,3] .= 0

    # Enforce free-slip no-flux boundary condition on velocity
    dUdt[bb,2], dUdt[bb,3] = 
        removeNormalComponent!(dUdt[bb,2], dUdt[bb,3], Tx, Ty)

    # # Dirichlet boundary condition for energy?  Still not sure.
    # dUdt[bb,4] .= 0

    return dUdt
    
end

