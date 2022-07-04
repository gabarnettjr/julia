
function rk3!(t::Float64, U::q, q1::q, q2::q, Utmp::q,
              pPrime::Array{Float64,1}, rhoPrime::Array{Float64,1},
              Tprime::Array{Float64,1}, D::Array{Float64,1},
              tmp1::Array{Float64,1}, tmp2::Array{Float64,1},
              ci::constantInputs)

    q1 = ODEfunction!(t, U, q1,
                      pPrime, rhoPrime, Tprime, D,
                      tmp1, tmp2, ci)

    t = t + ci.dt / 3.
    Utmp = qMult(ci.dt / 3., q1, Utmp)
    Utmp = qAdd(U, Utmp, Utmp)
    q2 = ODEfunction!(t, Utmp, q2,
                      pPrime, rhoPrime, Tprime, D,
                      tmp1, tmp2, ci)

    t = t + ci.dt / 3.
    Utmp = qMult(2. * ci.dt / 3., q2, Utmp)
    Utmp = qAdd(U, Utmp, Utmp)
    q2 = ODEfunction!(t, Utmp, q2,
                      pPrime, rhoPrime, Tprime, D,
                      tmp1, tmp2, ci)

    Utmp = qMult(3., q2, Utmp)
    Utmp = qAdd(q1, Utmp, Utmp)
    Utmp = qMult(ci.dt / 4., Utmp, Utmp)
    U = qAdd(U, Utmp, U)

    return  U

end

