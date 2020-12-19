
function ab3!(t, U, odefun!, dt, f1, f2, f3)

    U = U .+ dt ./ 12 .* (23 .* f3 .- 16 .* f2 .+ 5 .* f1)

    f1 = f2
    f2 = f3
    f3 = odefun!(t+dt, U, f3)

    return U, f1, f2, f3

end

