
function getInitialConditions(x, Cv, R, s, radius, dr)

    T_0 = 300 .* ones(size(x))
    e_0 = Cv .* T_0

    p_0 = 10^5 .* ones(size(x))
    
    q_0 = zeros(size(x))
    for i in 1 : length(x)
        r = sqrt((x[i] + 3/4*radius) ^ 2 + (y[i] - 0)^2)
        if r < 1
            q_0[i] = 1
        end
    end

    rho_0 = p_0 ./ R ./ T_0

    return rho_0, e_0, q_0, p_0

end
