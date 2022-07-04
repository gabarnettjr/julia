
function rk2!(t, U, odefun!, dt, q1)

    q1 = odefun!(t,       U,         q1)
    q1 = odefun!(t+dt./2, U+dt/2*q1, q1)

    U = U + dt * q1

    return U

end
