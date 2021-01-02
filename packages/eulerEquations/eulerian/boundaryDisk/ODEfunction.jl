
# include("removeNormalComponent.jl")
# include("fanPrime.jl")

function ODEfunction!(t, U, dUdt,
                      Wx, Wy, aWhv, rho_0, e_0, p_0,
                      Tx, Ty, bb, ff,
                      Cv, R)

    # Enforce free-slip no-flux boundary condition on velocity
    U[bb,2], U[bb,3] = removeNormalComponent!(U[bb,2], U[bb,3], Tx, Ty)

    # Equation of state for an ideal gas, p = rho * R * T
    p = U[:,1] .* R .* (U[:,4] ./ Cv)

    # Divergence of the velocity field (used twice)
    div = Wx * U[:,2] .+ Wy * U[:,3]

    # (d/dt)(rho)
    rho_0 = U[:,1] .- rho_0
    dUdt[:,1] = -U[:,2] .* (Wx * rho_0) .-
                 U[:,3] .* (Wy * rho_0) .-
                 U[:,1] .* div .+
                 aWhv * rho_0

    # (d/dt)(u)
    p_0 = p .- p_0
    dUdt[:,2] = -U[:,2] .* (Wx * U[:,2]) .- U[:,3] .* (Wy * U[:,2]) .-
                 (Wx * p_0) ./ U[:,1] .+
                 aWhv * U[:,2]

    # (d/dt)(v)
    dUdt[:,3] = -U[:,2] .* (Wx * U[:,3]) .- U[:,3] .* (Wy * U[:,3]) .-
                 (Wy * p_0) ./ U[:,1] .+
                 aWhv * U[:,3]
    
    # (d/dt)(e)
    e_0 = U[:,4] .- e_0
    dUdt[:,4] = -U[:,2] .* (Wx * e_0) .-
                 U[:,3] .* (Wy * e_0) .-
                 p .* div ./ U[:,1] .+
                 aWhv * e_0

    # Internal forcing (a fan is blowing air into the round room)
    dUdt[ff,2] .= fanPrime(t)
    dUdt[ff,3] .= 0

    # Enforce free-slip no-flux boundary condition on velocity (AGAIN!)
    dUdt[bb,2], dUdt[bb,3] = removeNormalComponent!(dUdt[bb,2], dUdt[bb,3],
                                                    Tx, Ty)

    # Dirichlet boundary condition for energy
    dUdt[bb,4] .= 0

    return dUdt
    
end

