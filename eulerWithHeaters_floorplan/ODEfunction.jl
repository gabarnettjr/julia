
# include("../packages/vectors/removeNormalComponent.jl")
# include("../packages/eulerEquations/fanPrime.jl")
# include("../packages/eulerEquations/heatPrime.jl")

function ODEfunction!(t::Float64, U::q, dUdt::q,
                      pPrime::Array{Float64,1}, rhoPrime::Array{Float64,1},
                      Tprime::Array{Float64,1}, D::Array{Float64,1},
                      tmp1::Array{Float64,1}, tmp2::Array{Float64,1},
                      ci::constantInputs)

    # use equation of state, then subtract background pressure
    pPrime .= U.rho .* ci.R .* U.T .- ci.p_0

    # Start out setting D equal to du/dx (updated later)
    D .= phs2_mult(ci.wx, U.u, ci.jj, D)

    # (d/dt)(u)
    dUdt.u .= .-U.u .* D .-
                U.v .* phs2_mult(ci.wy, U.u, ci.jj, Tprime) .-
                phs2_mult(ci.wx, pPrime, ci.jj, tmp1) ./ U.rho .+
                phs2_mult(ci.awhv, U.u, ci.jj, tmp2)

    # Temporarily set Tprime equal to dv/dy (half of D below)
    Tprime .= phs2_mult(ci.wy, U.v, ci.jj, Tprime)

    # (d/dt)(v)
    dUdt.v .= .-U.u .* phs2_mult(ci.wx, U.v, ci.jj, rhoPrime) .-
                U.v .* Tprime .-
                phs2_mult(ci.wy, pPrime, ci.jj, tmp1) ./ U.rho .+
                phs2_mult(ci.awhv, U.v, ci.jj, tmp2)

    # divergence of the velocity field
    D .= D .+ Tprime

    # (d/dt)(rho)
    rhoPrime .= U.rho .- ci.rho_0
    dUdt.rho .= .-U.u .* phs2_mult(ci.wx, rhoPrime, ci.jj, pPrime) .-
                  U.v .* phs2_mult(ci.wy, rhoPrime, ci.jj, tmp1) .-
                  U.rho .* D .+
                  phs2_mult(ci.awhv, rhoPrime, ci.jj, tmp2)
   
    # (d/dt)(T)
    Tprime .= U.T .- ci.T_0
    dUdt.T .= .-U.u .* phs2_mult(ci.wx, Tprime, ci.jj, pPrime) .-
                U.v .* phs2_mult(ci.wy, Tprime, ci.jj, tmp1) .-
                (ci.R ./ ci.Cv) .* U.T .* D .+
                phs2_mult(ci.awhv, Tprime, ci.jj, tmp2)
    
    # Fans blowing west to east
    dUdt.u[ci.ff2] .= fanPrime(t; fs = ci.fanSpeed)
    dUdt.v[ci.ff2] .= 0.
    # Fans blowing south to north
    dUdt.u[ci.ff3] .= 0.
    dUdt.v[ci.ff3] .= fanPrime(t; fs = ci.fanSpeed)
    # Fans blowing east to west
    dUdt.u[ci.ff4] .= -fanPrime(t; fs = ci.fanSpeed)
    dUdt.v[ci.ff4] .= 0.
    # Fans blowing north to south
    dUdt.u[ci.ff5] .= 0.
    dUdt.v[ci.ff5] .= -fanPrime(t; fs = ci.fanSpeed)

    # Internal forcing of heat (the fan is also a heater)
    dUdt.T[ci.ff2] .= heatPrime(t; td = tempDiff)
    dUdt.T[ci.ff3] .= heatPrime(t; td = tempDiff)
    dUdt.T[ci.ff4] .= heatPrime(t; td = tempDiff)
    dUdt.T[ci.ff5] .= heatPrime(t; td = tempDiff)

    # Enforce Dirichlet BC for the temperature (not sure about this)
    dUdt.T[ci.bb] .= 0.

    # # Enforce Neumann BC for the temperature (not sure about this)
    # dUdt.T[bb] = dUdt.T[bc]

    # Enforce free-slip no-flux boundary condition for velocity
    for i = 1 : length(bb)
        tmp = dUdt.u[bb[i]] * Tx[i] + dUdt.v[bb[i]] * Ty[i]
        dUdt.u[bb[i]] = tmp * Tx[i]
        dUdt.v[bb[i]] = tmp * Ty[i]
    end

    return dUdt
    
end

