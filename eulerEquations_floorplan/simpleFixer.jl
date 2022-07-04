
# The mass fixer works by estimating the average value of the
# density, and correcting it to the original average from t = 0.  It does the
# same thing for the tracer density.

function simpleFixer!(U, rho_0_avg, rhoq_0_avg, E_0_avg,
                      cs1, cs2, cs3, cs4, numSquares,
                      Cv, R, p_0)

    # Mass fixer
    rho = p_0 .* U[:,1] .^ (Cv/R) ./ R ./ U[:,4]
    rhoAvg = ((1/4)*sum(rho[cs1]) + (2/4)*sum(rho[cs2]) +
              (3/4)*sum(rho[cs3]) + (4/4)*sum(rho[cs4])) / numSquares
    rho = rho .+ (rho_0_avg .- rhoAvg)
    U[:,1] = (rho .* R .* U[:,4] ./ p_0) .^ (R/Cv)

    # Tracer mass fixer (and monotonicity fixer)
    q = U[:,5]
    for i in 1 : 4
        q[q .> 1] .= 1
        q[q .< 0] .= 0
        rhoq = rho .* q
        rhoqAvg = ((1/4)*sum(rhoq[cs1]) + (2/4)*sum(rhoq[cs2]) +
                   (3/4)*sum(rhoq[cs3]) + (4/4)*sum(rhoq[cs4])) / numSquares
        rhoq = rhoq .+ (rhoq_0_avg .- rhoqAvg)
        q = rhoq ./ rho
    end
    U[:,5] = q

    # # Total energy fixer (Can't use this anymore!)
    # uSquared = U[:,2] .^ 2 .+ U[:,3] .^ 2
    # e = Cv .* U[:,1] .* U[:,4]
    # E = rho .* (e .+ (1/2) .* uSquared)
    # Eavg = ((1/4)*sum(E[cs1]) + (2/4)*sum(E[cs2]) +
    #         (3/4)*sum(E[cs3]) + (4/4)*sum(E[cs4])) / numSquares
    # E = E .+ (E_0_avg .- Eavg)
    # e = E ./ rho - (1/2) .* uSquared
    # U[:,4] = e ./ Cv ./ U[:,1]

    return U

end

