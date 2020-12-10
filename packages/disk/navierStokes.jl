###############################################################################

function getConstants()
    
    mu = 1.8 * 10^-5
    # mu = 0
    
    k = 0
    lam = k - 2/3 * mu
    Cv = 717
    R = 287
    
    return mu, k, lam, Cv, R
    
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

function ODEfunction!(t, U, dUdt,
x, y, Wx, Wy, aWhv, rho_0, e_0, p_0, bb1, bb2, bb3, mu, k, lam, Cv, R)
    
    T = U[:,4] ./ Cv
    p = U[:,1] .* R .* T

    # Boundary conditions
    U[bb1, 2] = exp.(-9 .* y[bb1] .^ 2) .* 10 .* t.^3 ./ (20 .+ t.^3)
    U[bb1, 3] .= 0
    U[bb2, 2] = -y[bb2] ./ x[bb2] .* U[bb2, 3]
    U[bb3, 3] = -x[bb3] ./ y[bb3] .* U[bb3, 2]

    dudx = Wx * U[:,2]
    dudy = Wy * U[:,2]
    dvdx = Wx * U[:,3]
    dvdy = Wy * U[:,3]
    
    dUdt[:,1] = -U[:,2] .* (Wx * (U[:,1] - rho_0)) .-
                 U[:,3] .* (Wy * (U[:,1] - rho_0)) .-
                 U[:,1] .* (dudx .+ dvdy) .+
                 aWhv * (U[:,1] - rho_0)
    
    dUdt[:,2] = -U[:,2] .* dudx .- U[:,3] .* dudy .-
                 1 ./ U[:,1] .* (Wx * (p - p_0) .-
                 Wx * (lam .* (dudx .+ dvdy)) .-
                 Wx * (mu .* 2 .* dudx) .-
                 Wy * (mu .* (dvdx .+ dudy)))
                 # aWhv * U[:,2]
    
    dUdt[:,3] = -U[:,2] .* dvdx .- U[:,3] .* dvdy .-
                 1 ./ U[:,1] .* (Wy * (p - p_0) .-
                 Wy * (lam .* (dudx .+ dvdy)) .-
                 Wx * (mu .* (dudy .+ dvdx)) .-
                 Wy * (mu .* 2 .* dvdy))
                 # aWhv * U[:,3]
    
    dUdt[:,4] = -U[:,2] .* (Wx * (U[:,4] - e_0)) .-
                 U[:,3] .* (Wy * (U[:,4] - e_0)) .-
                 1 ./ U[:,1] .* (p .* (dudx .+ dvdy) .-
                 # Wx * (k .* (Wx*T)) .-
                 # Wy * (k .* (Wy*T)) .-
                 lam .* (dudx.^2 .+ dvdy.^2)) .+
                 mu ./ U[:,1] .* (2 .* dudx.^2 .+
                 (dudy .+ dvdx) .* dvdx .+
                 (dvdx .+ dudy) .* dudy .+
                 2 .* dvdy.^2) .+
                 aWhv * (U[:,4] - e_0)

    # in-flow boundary condition
    dUdt[bb1, 2] = exp.(-9 .* y[bb1] .^ 2) .* 600 .* t.^2 ./ (20 .+ t.^3) .^ 2
    dUdt[bb1, 3] .= 0
    # dUdt[bb1, 2] .= 0

    # no-flux boundary condition
    dUdt[bb2, 2] = -y[bb2] ./ x[bb2] .* dUdt[bb2, 3]
    dUdt[bb3, 3] = -x[bb3] ./ y[bb3] .* dUdt[bb3, 2]
    # dUdt[bb2, 2:3] .= 0
    # dUdt[bb3, 2:3] .= 0

    # Dirichlet boundary condition
    # dUdt[bb2, 1] .= 0
    # dUdt[bb3, 1] .= 0
    # dUdt[bb2, 4] .= 0
    # dUdt[bb3, 4] .= 0

    return dUdt
    
end

###############################################################################
