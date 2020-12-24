
function getInitialConditions(x, Cv, R, domain)

    T_0 = 300 .* ones(size(x))
    e_0 = Cv .* T_0

    if domain == 1

        p_0 = 10^5 .+ 100 .* (10 .- x)

    elseif domain == 2

        p_0 = 10^5 .* ones(size(x))
        for i = 1 : length(x)
            if (x[i] < 4)
                p_0[i] = 10^5 + 100 * (4 - x[i])
            end
        end

    end

    rho_0 = p_0 ./ R ./ T_0

    u_0 = zeros(size(x))
    v_0 = zeros(size(x))
    
    return rho_0, u_0, v_0, e_0, p_0
    
end

