###########################################################################

function rk1!(t, U, odefun, dt, q1)

    q1 = odefun(t, U, q1)

    t = t + dt
    U = U + dt .* q1

    return t, U

end

###########################################################################

function rk2!(t, U, odefun, dt, q1)

    q1 = odefun(t,       U,         q1)
    q1 = odefun(t+dt./2, U+dt/2*q1, q1)

    t = t + dt
    U = U + dt * q1

    return t, U

end

###########################################################################

function rk3!(t, U, odefun, dt, q1, q2)

    q1 = odefun(t,        U,           q1)
    q2 = odefun(t+dt/3,   U+dt/3*q1,   q2)
    q2 = odefun(t+2*dt/3, U+2*dt/3*q2, q2)

    t = t + dt
    U = U + dt/4 * (q1 + 3*q2)

    return t, U

end

###########################################################################

function rk4!(t, U, odefun, dt, q1, q2, q3, q4)

    q1 = odefun(t,      U,         q1)
    q2 = odefun(t+dt/2, U+dt/2*q1, q2)
    q3 = odefun(t+dt/2, U+dt/2*q2, q3)
    q4 = odefun(t+dt,   U+dt*q3,   q4)

    t = t + dt
    U = U + dt/6 * (q1 + 2*q2 + 2*q3 + q4)

    return t, U

end

###########################################################################
