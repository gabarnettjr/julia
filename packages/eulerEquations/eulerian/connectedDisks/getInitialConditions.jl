
function getInitialConditions(x, Cv, R, radius, dr)
    
    T_0 = 300 .* ones(length(x))
    
    e_0 = Cv .* T_0
    
    # p_0 = 10^5 .* ones(length(x))
    # for i in 1 : length(p_0)
    #     r = sqrt((x[i] + (2*radius+2*dr)) ^ 2 + y[i]^2)
    #     if (r < radius)
    #         p_0[i] = p_0[i] + 500
    #     end
    # end

    p_0 = 10^5 .- 100 .* x

    # r = sqrt.(x .^ 2 .+ y .^ 2)
    # p_0 = 10^5 .+ 100 .* r

    rho_0 = p_0 ./ R ./ T_0
    
    u_0 = zeros(length(x))
    
    v_0 = zeros(length(x))
    
    return rho_0, u_0, v_0, e_0, p_0
    
end

