
# include("freeSlipNoFlux.jl")

function ODEfunction!(t, U, dUdt,
                      Wx, Wy, aWhv, rho_0, e_0,
                      i_n1, j_n1, i_n2, j_n2, i_n3, j_n3, i_n4, j_n4,
                      i_n12, j_n12, i_n14, j_n14, i_n23, j_n23, i_n34, j_n34,
                      bool_all, nanMats, zeroVecs,
                      Cv, R)
    
    U = vecsToMats!(U, bool_all, nanMats)

    U = freeSlipNoFlux!(U, i_n1, j_n1, i_n2, j_n2, i_n3, j_n3, i_n4, j_n4,
                        i_n12, j_n12, i_n14, j_n14, i_n23, j_n23, i_n34, j_n34)

    U = matsToVecs!(U, bool_all, zeroVecs)

    p = U[:,1] .* R .* (U[:,4] ./ Cv)

    dudx = Wx * U[:,2]
    dvdy = Wy * U[:,3]
    
    dUdt[:,1] = -U[:,2] .* (Wx * U[:,1]) .-
                 U[:,3] .* (Wy * U[:,1]) .-
                 U[:,1] .* (dudx .+ dvdy) .+
                 aWhv * (U[:,1] - rho_0)
    
    dUdt[:,2] = -U[:,2] .* dudx .- U[:,3] .* (Wy * U[:,2]) .-
                 (Wx * p) ./ U[:,1] .+
                 aWhv * U[:,2]
    
    dUdt[:,3] = -U[:,2] .* (Wx * U[:,3]) .- U[:,3] .* dvdy .-
                 (Wy * p) ./ U[:,1] .+
                 aWhv * U[:,3]
    
    dUdt[:,4] = -U[:,2] .* (Wx * U[:,4]) .-
                 U[:,3] .* (Wy * U[:,4]) .-
                 p .* (dudx .+ dvdy) ./ U[:,1] .+
                 aWhv * (U[:,4] - e_0)

    return dUdt
    
end

