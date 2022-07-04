
# The initial conditions for energy and pressure are constant, and the initial
# condition for density is defined using the ideal gas law, p = rho * R * T.
# The initial value of the tracer q is zero everywhere, except some area,
# usually close to where a fan is blowing in the spatial domain.  The tracer
# is so you can watch what happens to a particular patch of air as time passes.

function getInitialConditions(R, bgTemp)

    T_0 = bgTemp .* ones(size(x))
    # e_0 = Cv .* T_0

    p_0 = 10^5 .* ones(size(x))
    
    rho_0 = p_0 ./ R ./ T_0

    return rho_0, p_0, T_0

end
