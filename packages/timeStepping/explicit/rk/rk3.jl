
function rk3!(t, U, odefun!, dt, q1, q2)

    q1 = odefun!(t,        U,           q1)
    q2 = odefun!(t+dt/3,   U+dt/3*q1,   q2)
    q2 = odefun!(t+2*dt/3, U+2*dt/3*q2, q2)

    U = U + dt/4 * (q1 + 3*q2)

    return U

end
