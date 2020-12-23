
function getInitialConditions(x, Cv, R)

    T_0 = 300 .* ones(size(x))
    e_0 = Cv .* T_0
    p_0 = 10^5 .- 100 .* x
    rho_0 = p_0 ./ R ./ T_0

    # p_0 = 10^5 .- 100 .* x
    # rho_0 = 1.16144 .* ones(size(x))
    # T_0 = p_0 ./ rho_0 ./ R
    # e_0 = Cv .* T_0

    u_0 = zeros(size(x))
    v_0 = zeros(size(x))
    
    return rho_0, u_0, v_0, e_0, p_0
    
end

