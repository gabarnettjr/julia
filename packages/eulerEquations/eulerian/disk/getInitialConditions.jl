
function getInitialConditions(x, Cv, R)
    
    # T_0 = 300 .* ones(length(x))
    
    # e_0 = Cv .* T_0
    
    # p_0 = 10^5 .* ones(length(x))

    # p_0 = 10^5 .- 100 .* x
    # rho_0 = 1.16144 .* ones(size(x))
    # T_0 = p_0 ./ rho_0 ./ R
    # e_0 = Cv .* T_0

    T_0 = 300 .* ones(size(x))
    e_0 = Cv .* T_0
    p_0 = 10^5 .- 100 .* x
    rho_0 = p_0 ./ R ./ T_0

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
    
    # rho_0 = p_0 ./ R ./ T_0
    
    u_0 = zeros(size(x))
    
    v_0 = zeros(size(x))
    
    return rho_0, u_0, v_0, e_0, p_0
    
end

