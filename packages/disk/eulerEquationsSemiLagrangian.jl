###############################################################################

function getConstants()
    
    Cv = 717
    R = 287
    
    return Cv, R
    
end

###############################################################################

function getInitialConditions(x, Cv, R)
    
    T_0 = 300 .* ones(length(x))
    
    e_0 = Cv .* T_0
    
    p_0 = 10^5 .* ones(length(x))
    
    rho_0 = p_0 ./ R ./ T_0
    
    u_0 = zeros(length(x))
    
    v_0 = zeros(length(x))
    
    return rho_0, u_0, v_0, e_0, p_0
    
end

###############################################################################

function appendGhostNodes!(x, y, dr)

    xb = []
    yb = []

    for i in 1 : length(x)
        if abs(sqrt(x[i]^2 + y[i]^2) - 1) < 1e-3
            xb = vcat(xb, x[i])
            yb = vcat(yb, y[i])
        end
    end

    th = atan.(yb, xb)

    x = vcat(x, xb .+ dr .* cos.(th))
    y = vcat(y, yb .+ dr .* sin.(th))

    x = vcat(x, xb .+ 2 .* dr .* cos.(th))
    y = vcat(y, yb .+ 2 .* dr .* sin.(th))

    return x, y

end

###############################################################################

function getIndices(x, y, dr, radius)
    
    ind_ghost = []
    ind_noGhost = []
    ind_outer = []
    ind_fan = []

    nLayers = Int(round((radius + dr/2) / dr / 2))
    fanRad = radius + dr/2 - nLayers * dr
    
    for i in 1 : length(x)
    
        r = sqrt(x[i] ^2 + y[i] ^2)
        th = atan(y[i], x[i])
        
        if abs(r - (radius+dr/2)) < 1e-3
            ind_ghost = vcat(ind_ghost, i)
        elseif abs(r - (radius-dr/2)) < 1e-3
            ind_outer = vcat(ind_outer, i)
        elseif abs(r - fanRad) < 1e-3 && (th >= 11*pi/12 || th <= -11*pi/12)
            ind_fan = vcat(ind_fan, i)
        end

        if r < radius
            ind_noGhost = vcat(ind_noGhost, i)
        end
    
    end

    return ind_ghost, ind_noGhost, ind_outer, ind_fan

end

###############################################################################

function interp!(U, x, y, phs, pol, stc)

    W = getDM(hcat(x, y)', U[:,5:6]', [0 0], phs, pol, stc, 0)[1]
    U[:,1:4] = W * U[:,1:4]
    U[:,5] = x
    U[:,6] = y

    return U

end

###############################################################################

function freeSlipNoFlux!(U, ind_ghost, ind_noGhost, ind_outer, radius,
phs, pol, stc, e_0)

    # Modify the location of the ghost nodes so they are reflected over the
    # boundary from the outer nodes:
    dr = radius .- sqrt.(U[ind_outer,5] .^2 .+ U[ind_outer,6] .^2)
    th = atan.(U[ind_outer,6], U[ind_outer,5])
    U[ind_ghost,5] = U[ind_outer,5] .+ 2 .* dr .* cos.(th)
    U[ind_ghost,6] = U[ind_outer,6] .+ 2 .* dr .* sin.(th)

    # Get the extrapolation weights to go from outer nodes to ghost nodes:
    w, jj = getDM(U[ind_ghost,5:6]', U[ind_noGhost,5:6]',
                  [0 0], phs, pol, stc, 0; getSparse = false)[[1,2]]

    # Extrapolate density to the ghost nodes:
    U[ind_ghost,1] = sum(w .* U[jj,1], dims = 1)

    # Enforce Dirichlet boundary condition on the energy:
    U[ind_ghost,4] = -U[ind_outer,4] .+ 2 .* e_0[ind_outer]

    # repeat theta at each outer node so it can be multiplied by stencils:
    thth = repeat(th', stc)

    # Get the tangent velocity on the stencils near the boundary:
    uth_outer = -sin.(thth) .* U[jj,2] .+ cos.(thth) .* U[jj,3]

    # Extrapolation of the tangent velocity to the ghost nodes:
    uth_ghost = sum(w .* uth_outer, dims = 1)'
    
    # Get the normal velocity on the outer-most non-ghost level:
    ur_outer = cos.(th) .* U[ind_outer,2] .+ sin.(th) .* U[ind_outer,3]

    # Set ghost nodes values to enforce no-flux condition:
    ur_ghost = -ur_outer

    # Set the values of u and v on the ghost nodes:
    U[ind_ghost,2] = cos.(th) .* ur_ghost .- sin.(th) .* uth_ghost
    U[ind_ghost,3] = sin.(th) .* ur_ghost .+ cos.(th) .* uth_ghost

    return U

end

###############################################################################

function fan!(t, U, ind_fan, phs, pol, stc)
    
    f = 10 .* t .^ 3 ./ (1 .+ t .^ 3)
    U[ind_fan,2] .= f
    U[ind_fan,3] .= 0

    W = getDM(U[ind_fan,5:6]', U[:,5:6]', [0 0], phs, pol, stc, 0)[1]
    U[ind_fan,[1,4]] = W * U[:,[1,4]]

    return U

end

###############################################################################

function ODEfunction!(t, U, dUdt,
radius, rho_0, e_0, p_0, phs, pol, stc, K,
ind_ghost, ind_noGhost, ind_outer, ind_fan,
Cv, R, a)
    
    U = freeSlipNoFlux!(U, ind_ghost, ind_noGhost, ind_outer, radius,
                        phs, pol, stc, e_0)
    # U = fan!(t, U, ind_fan, phs, pol, stc)

    # Get the differentiation matrices on the (possibly) moved nodes:
    Wx, tr, ind_nn =
        getDM(U[:,5:6]', U[:,5:6]', [1 0], phs, pol, stc, 0)[[1,3,4]]
    Wy = getDM(U[:,5:6]', U[:,5:6]', [0 1], phs, pol, stc, 0;
               tree = tr, idx = ind_nn)[1]
    aWhv = a .* getDM(U[:,5:6]', U[:,5:6]', [-1 -1], phs, pol, stc, K;
                      tree = tr, idx = ind_nn)[1]

    p = U[:,1] .* R .* (U[:,4] ./ Cv)

    div = Wx * U[:,2] .+ Wy * U[:,3]
    
    dUdt[:,1] = -U[:,1] .* div .+
                aWhv * (U[:,1] - rho_0)
    
    dUdt[:,2] = -(Wx * (p - p_0)) ./ U[:,1] .+
                aWhv * U[:,2]
    
    dUdt[:,3] = -(Wy * (p - p_0)) ./ U[:,1] .+
                aWhv * U[:,3]
    
    dUdt[:,4] = -p .* div ./ U[:,1] .+
                aWhv * (U[:,4] - e_0)

    dUdt[:,5] = U[:,2]

    dUdt[:,6] = U[:,3]

    # dUdt[ind_fan,:] .= 0

    return dUdt
    
end

###############################################################################
