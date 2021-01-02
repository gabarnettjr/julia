
function simpleFixer!(U, q, rho_0_avg, rhoq_0_avg, e_0_avg, ii)

    n = size(U, 1)
    
    # Mass fixer
    rhoAvg = sum(U[:,1]) / n
    U[:,1] = U[:,1] .+ (rho_0_avg .- rhoAvg)

    # Tracer mass fixer (and monotonicity fixer)
    for i in 1 : 4
        q[q .> 1] .= 1
        q[q .< 0] .= 0
        rhoq = U[:,1] .* q
        rhoqAvg = sum(rhoq) / n
        rhoq = rhoq .+ (rhoq_0_avg .- rhoqAvg)
        q = rhoq ./ U[:,1]
    end
    q[q .> 1] .= 1
    q[q .< 0] .= 0

    # Horizontal momentum fixer
    rhoU = U[ii,1] .* U[ii,2]
    rhoUavg = sum(rhoU) / length(rhoU)
    rhoU = rhoU .- rhoUavg
    U[ii,2] = rhoU ./ U[ii,1]

    # Vertical momentum fixer
    rhoV = U[ii,1] .* U[ii,3]
    rhoVavg = sum(rhoV) / length(rhoV)
    rhoV = rhoV .- rhoVavg
    U[ii,3] = rhoV ./ U[ii,1]

    # Total energy fixer
    E_0_avg = rho_0_avg * e_0_avg
    uSquared = U[:,2] .^ 2 .+ U[:,3] .^ 2
    E = U[:,1] .* (U[:,4] .+ (1/2) .* uSquared)
    Eavg = sum(E) / n
    E = E .+ (E_0_avg .- Eavg)
    U[:,4] = E ./ U[:,1] - (1/2) .* uSquared

    return U, q

end

