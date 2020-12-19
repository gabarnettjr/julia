
function rk4!(t, U, odefun!, dt, q1, q2, q3, q4)

    q1 = odefun!(t,      U,         q1)
    q2 = odefun!(t+dt/2, U+dt/2*q1, q2)
    q3 = odefun!(t+dt/2, U+dt/2*q2, q3)
    q4 = odefun!(t+dt,   U+dt*q3,   q4)

    U = U + dt/6 * (q1 + 2*q2 + 2*q3 + q4)

    return U

end
