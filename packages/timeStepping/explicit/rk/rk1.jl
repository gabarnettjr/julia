
function rk1!(t, U, odefun!, dt, q1)

    q1 = odefun!(t, U, q1)

    U = U + dt .* q1

    return U

end
