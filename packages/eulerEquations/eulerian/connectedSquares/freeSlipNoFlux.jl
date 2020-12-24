
function freeSlipNoFlux!(U, i_n1, j_n1, i_n2, j_n2, i_n3, j_n3, i_n4, j_n4,
i_n12, j_n12, i_n14, j_n14, i_n23, j_n23, i_n34, j_n34)

    for k in 1 : length(i_n1)
        # extrapolate rho, v, and e:
        U[i_n1[k], j_n1[k], [1,3,4]] = 2 .* U[i_n1[k], j_n1[k] + 1, [1,3,4]] .-
                                            U[i_n1[k], j_n1[k] + 2, [1,3,4]]
        # enforce no flux on u:
        U[i_n1[k], j_n1[k], 2] = -U[i_n1[k], j_n1[k] + 1, 2]
    end

    for k in 1 : length(i_n2)
        # extrapolate rho, u, and e:
        U[i_n2[k], j_n2[k], [1,2,4]] = 2 .* U[i_n2[k] + 1, j_n2[k], [1,2,4]] .-
                                            U[i_n2[k] + 2, j_n2[k], [1,2,4]]
        # enforce no flux on v:
        U[i_n2[k], j_n2[k], 3] = -U[i_n2[k] + 1, j_n2[k], 3]
    end

    for k in 1 : length(i_n3)
        # extrapolate rho, v, and e:
        U[i_n3[k], j_n3[k], [1,3,4]] = 2 .* U[i_n3[k], j_n3[k] - 1, [1,3,4]] .-
                                            U[i_n3[k], j_n3[k] - 2, [1,3,4]]
        # enforce no flux on u:
        U[i_n3[k], j_n3[k], 2] = -U[i_n3[k], j_n3[k] - 1, 2]
    end

    for k in 1 : length(i_n4)
        # extrapolate rho, u, and e:
        U[i_n4[k], j_n4[k], [1,2,4]] = 2 .* U[i_n4[k] - 1, j_n4[k], [1,2,4]] .-
                                            U[i_n4[k] - 2, j_n4[k], [1,2,4]]
        # enforce no flux on v:
        U[i_n4[k], j_n4[k], 3] = -U[i_n4[k] - 1, j_n4[k], 3]
    end

    for k in 1 : length(i_n12)
        # # extrapolate rho and e from both sides and average:
        # tmp1 = 2 .* U[i_n12[k] + 1, j_n12[k],     [1,4]] .-
        #             U[i_n12[k] + 2, j_n12[k],     [1,4]]
        # tmp2 = 2 .* U[i_n12[k],     j_n12[k] + 1, [1,4]] .-
        #             U[i_n12[k],     j_n12[k] + 2, [1,4]]
        # U[i_n12[k], j_n12[k], [1,4]] = (tmp1 .+ tmp2) ./ 2
        # just extrapolate diagonally:
        U[i_n12[k], j_n12[k], [1,4]] =
            2 .* U[i_n12[k] + 1, j_n12[k] + 1, [1,4]] .-
                 U[i_n12[k] + 2, j_n12[k] + 2, [1,4]]
        # enforce no flux on u:
        U[i_n12[k], j_n12[k], 2] = -U[i_n12[k],     j_n12[k] + 1, 2]
        # enforce no flux on v:
        U[i_n12[k], j_n12[k], 3] = -U[i_n12[k] + 1, j_n12[k],     3]
    end

    for k in 1 : length(i_n14)
        # extrapolate rho and e:
        # tmp1 = 2 .* U[i_n14[k] - 1, j_n14[k],     [1,4]] .-
        #             U[i_n14[k] - 2, j_n14[k],     [1,4]]
        # tmp2 = 2 .* U[i_n14[k],     j_n14[k] + 1, [1,4]] .-
        #             U[i_n14[k],     j_n14[k] + 2, [1,4]]
        # U[i_n14[k], j_n14[k], [1,4]] = (tmp1 .+ tmp2) ./ 2
        U[i_n14[k], j_n14[k], [1,4]] =
            2 .* U[i_n14[k] - 1, j_n14[k] + 1, [1,4]] .-
                 U[i_n14[k] - 2, j_n14[k] + 2, [1,4]]
        # enforce no flux on u:
        U[i_n14[k], j_n14[k], 2] = -U[i_n14[k],     j_n14[k] + 1, 2]
        # enforce no flux on v:
        U[i_n14[k], j_n14[k], 3] = -U[i_n14[k] - 1, j_n14[k],     3]
    end

    for k in 1 : length(i_n23)
        # extrapolate rho and e:
        # tmp1 = 2 .* U[i_n23[k] + 1, j_n23[k],     [1,4]] .-
        #             U[i_n23[k] + 2, j_n23[k],     [1,4]]
        # tmp2 = 2 .* U[i_n23[k],     j_n23[k] - 1, [1,4]] .-
        #             U[i_n23[k],     j_n23[k] - 2, [1,4]]
        # U[i_n23[k], j_n23[k], [1,4]] = (tmp1 .+ tmp2) ./ 2
        U[i_n23[k], j_n23[k], [1,4]] =
            2 .* U[i_n23[k] + 1, j_n23[k] - 1, [1,4]] .-
                 U[i_n23[k] + 2, j_n23[k] - 2, [1,4]]
        # enforce no flux on u:
        U[i_n23[k], j_n23[k], 2] = -U[i_n23[k],     j_n23[k] - 1, 2]
        # enforce no flux on v:
        U[i_n23[k], j_n23[k], 3] = -U[i_n23[k] + 1, j_n23[k],     3]
    end

    for k in 1 : length(i_n34)
        # extrapolate rho and e:
        # tmp1 = 2 .* U[i_n34[k] - 1, j_n34[k],     [1,4]] .-
        #             U[i_n34[k] - 2, j_n34[k],     [1,4]]
        # tmp2 = 2 .* U[i_n34[k],     j_n34[k] - 1, [1,4]] .-
        #             U[i_n34[k],     j_n34[k] - 2, [1,4]]
        # U[i_n34[k], j_n34[k], [1,4]] = (tmp1 .+ tmp2) ./ 2
        U[i_n34[k], j_n34[k], [1,4]] =
            2 .* U[i_n34[k] - 1, j_n34[k] - 1, [1,4]] .-
                 U[i_n34[k] - 2, j_n34[k] - 2, [1,4]]
        # enforce no flux on u:
        u = -U[i_n34[k],     j_n34[k] - 1, 2]
        # enforce no flux on v:
        v = -U[i_n34[k] - 1, j_n34[k],     3]
        # enforce diagonal flux being zero (this might be more important)
        U[i_n34[k], j_n34[k], 2] = u - 1/2 * (u + v)
        U[i_n34[k], j_n34[k], 3] = v - 1/2 * (u + v)
    end

    return U

end

