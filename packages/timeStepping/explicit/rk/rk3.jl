
function rk3!(t, U, odefun!, dt, q1, q2, Utmp)

    q1 = odefun!(t, U, q1)

    Utmp = U + dt/3 * q1
    q2 = odefun!(t + dt/3, Utmp, q2)

    Utmp = U + 2*dt/3 * q2
    q2 = odefun!(t + 2*dt/3, Utmp, q2)

    Utmp = U + dt/4 * (q1 + 3 * q2)

    return Utmp

end
