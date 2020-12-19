###############################################################################

function getConstants()
    
    mu = 1.8 * 10^-5
    
    k = 0
    lam = k - 2/3 * mu
    Cv = 717
    R = 287
    
    return mu, k, lam, Cv, R
    
end

###############################################################################

function getInitialConditions(x, Cv, R)
    
    T_0 = 300 .* ones(length(x))
    # T_0 = 15 .* (20 .- x)
    
    e_0 = Cv .* T_0
    
    # p_0 = 10^5 .* ones(length(x))

    p_0 = 10^5 .- 100 .* x

    # r = sqrt.(x .^ 2 .+ y .^ 2)
    # p_0 = 10^5 .+ 100 .* r

    # p_0 = zeros(length(x))
    # for i = 1 : length(x)
    #     if x[i] < 0
    #         p_0[i] = 10^5 + 100
    #     elseif abs(x[i]) < 1e-3
    #         p_0[i] = 10^5
    #     else
    #         p_0[i] = 10^5 - 100
    #     end
    # end
    
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

function getIndices(x, y, dr, radius, thetaIn, thetaOut, outAngle)
    
    indA = []
    indB = []
    indC = []
    indA_inflow = []
    indB_inflow = []
    indC_inflow = []
    indA_outflow = []
    indB_outflow = []
    indC_outflow = []
    indFan = []
    
    nLayers = Int(round((radius + dr/2) / dr / 4))
    fanRad = radius + dr/2 - nLayers * dr
    
    for i in 1 : length(x)
    
        r = sqrt(x[i] ^2 + y[i] ^2)
        th = atan(y[i], x[i])
        
        if abs(r - (radius+dr/2)) < 1e-3

            if (th > pi - thetaIn/2 || th < -pi + thetaIn/2)
                indC_inflow = vcat(indC_inflow, i)
            elseif (th >= outAngle + thetaOut/2 || th <= outAngle - thetaOut/2)
                indC = vcat(indC, i)
            else
                indC_outflow = vcat(indC_outflow, i)
            end
            
        elseif abs(r - (radius-dr/2)) < 1e-3

            if (th > pi - thetaIn/2 || th < -pi + thetaIn/2)
                indB_inflow = vcat(indB_inflow, i)
            elseif (th >= outAngle + thetaOut/2 || th <= outAngle - thetaOut/2)
                indB = vcat(indB, i)
            else
                indB_outflow = vcat(indB_outflow, i)
            end

        elseif abs(r - (radius-3*dr/2)) < 1e-3

            if (th > pi - thetaIn/2 || th < -pi + thetaIn/2)
                indA_inflow = vcat(indA_inflow, i)
            elseif (th >= outAngle + thetaOut/2 || th <= outAngle - thetaOut/2)
                indA = vcat(indA, i)
            else
                indA_outflow = vcat(indA_outflow, i)
            end
        
        elseif abs(r - fanRad) < 1e-3 && (th >= 11*pi/12 || th <= -11*pi/12)

            indFan = vcat(indFan, i)

        end
    
    end

    return indA,         indB,         indC,
           indA_inflow,  indB_inflow,  indC_inflow,
           indA_outflow, indB_outflow, indC_outflow,
           indFan

end

###############################################################################

function freeSlipNoFlux!(U, th, indA, indB, indC, e_0)

    # Start by just extrapolating the density and energy to the ghost nodes:
    U[indC,[1,4]] = 2 .* U[indB,[1,4]] .- U[indA,[1,4]]

    # # Enforce a Dirichlet boundary condition on the energy:
    # U[indC,4] = 2 .* e_0[indC] - U[indB,4]

    # Get the tangent velocity on the two non-ghost levels:
    uthA = -sin.(th) .* U[indA,2] .+ cos.(th) .* U[indA,3]
    uthB = -sin.(th) .* U[indB,2] .+ cos.(th) .* U[indB,3]

    # Linear extrapolation of the tangent velocity to the ghost nodes:
    uthC = 2 .* uthB .- uthA
    
    # Get the normal velocity on the outer-most non-ghost level:
    urB = cos.(th) .* U[indB,2] .+ sin.(th) .* U[indB,3]

    # Set ghost nodes values to enforce no-flux condition:
    urC = -urB
    
    # Set the values of u and v on the ghost nodes:
    U[indC,2] = cos.(th) .* urC .- sin.(th) .* uthC
    U[indC,3] = sin.(th) .* urC .+ cos.(th) .* uthC

    return U

end

###############################################################################

function inflow!(t, U, y_inflow, ymax, indA_inflow, indB_inflow, indC_inflow)

    U[indC_inflow,:] = 2 .* U[indB_inflow,:] .- U[indA_inflow,:]
    
    f = 10 .* t .^ 3 ./ (1 .+ t .^ 3)

    # f = (exp.(-1 .* y_inflow .^ 2) .- exp(-1 .* ymax .^ 2)) .* f

    U[indC_inflow,2] .= 2 .* f .- U[indB_inflow,2]

    U[indC_inflow,3] .= -U[indB_inflow,3]

    return U

end

###############################################################################

function outflow!(U, indA_outflow, indB_outflow, indC_outflow)

    U[indC_outflow,:] = 2 .* U[indB_outflow,:] .- U[indA_outflow,:]

    return U

end

###############################################################################

function fan!(t, U, y_inflow, ymax, indFan)
    
    f = 10 .* t .^ 3 ./ (1 .+ t .^ 3)
    # f = (exp.(-1 .* y_inflow .^ 2) .- exp(-1 .* ymax .^ 2)) .* f

    U[indFan,2] .= f
    U[indFan,3] .= 0

    return U

end

###############################################################################

function ODEfunction!(t, U, dUdt,
x, y, th, y_inflow, ymax, radius, Wx, Wy, aWhv, rho_0, e_0, p_0,
indA, indB, indC,
indA_inflow, indB_inflow, indC_inflow,
indA_outflow, indB_outflow, indC_outflow,
indFan,
mu, k, lam, Cv, R)
    
    U = freeSlipNoFlux!(U, th, indA, indB, indC, e_0)
    # U = fan!(t, U, y_inflow, ymax, indFan)
    # U = inflow!(t, U, y_inflow, ymax, indA_inflow, indB_inflow, indC_inflow)
    # U = outflow!(U, indA_outflow, indB_outflow, indC_outflow)

    T = U[:,4] ./ Cv
    p = U[:,1] .* R .* T

    dudx = Wx * U[:,2]
    dudy = Wy * U[:,2]
    dvdx = Wx * U[:,3]
    dvdy = Wy * U[:,3]
    
    dUdt[:,1] = -U[:,2] .* (Wx * U[:,1]) .-
                 U[:,3] .* (Wy * U[:,1]) .-
                 U[:,1] .* (dudx .+ dvdy) .+
                 aWhv * (U[:,1] - rho_0)
    
    dUdt[:,2] = -U[:,2] .* dudx .- U[:,3] .* dudy .-
                 (Wx * p) ./ U[:,1] .+
                 # (Wx * (lam .* (dudx .+ dvdy)) .+
                 # Wx * (mu .* 2 .* dudx) .+
                 # Wy * (mu .* (dvdx .+ dudy))) ./ U[:,1] .+
                 aWhv * U[:,2]
    
    dUdt[:,3] = -U[:,2] .* dvdx .- U[:,3] .* dvdy .-
                 (Wy * p) ./ U[:,1] .+
                 # (Wy * (lam .* (dudx .+ dvdy)) .+
                 # Wx * (mu .* (dudy .+ dvdx)) .+
                 # Wy * (mu .* 2 .* dvdy)) ./ U[:,1] .+
                 aWhv * U[:,3]
    
    dUdt[:,4] = -U[:,2] .* (Wx * U[:,4]) .-
                 U[:,3] .* (Wy * U[:,4]) .-
                 p .* (dudx .+ dvdy) ./ U[:,1] .+
                 # (Wx * (k .* (Wx*T)) .+
                 # Wy * (k .* (Wy*T)) .+
                 # lam .* (dudx.^2 .+ dvdy.^2)) ./ U[:,1] .+
                 # mu .* (2 .* dudx.^2 .+
                 # (dudy .+ dvdx) .* dvdx .+
                 # (dvdx .+ dudy) .* dudy .+
                 # 2 .* dvdy.^2) ./ U[:,1] .+
                 aWhv * (U[:,4] - e_0)

    return dUdt
    
end

###############################################################################
