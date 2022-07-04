
function expandFloorplanMatrix(M::Array{Int64,2}, k::Int64)

    (m, n) = size(M)

    Mexp = zeros(Int, k*m, k*n)

    for j in 1 : n
        for i in 1 : m
            ii = k*(i - 1) + 1 : k*i
            jj = k*(j - 1) + 1 : k*j
            if i == 1 || i == m || j == 1 || j == n
                if i == 1
                    if j != 1 && j != n
                        tmp = ones(k,k)
                        tmp[1, :] .= 0
                        Mexp[ii, jj] = tmp
                    elseif j == 1
                        tmp = ones(k,k)
                        tmp[:, 1] .= 0
                        tmp[1, :] .= 0
                        Mexp[ii, jj] = tmp
                    elseif j == n
                        tmp = ones(k,k)
                        tmp[1, :] .= 0
                        tmp[:, k] .= 0
                        Mexp[ii, jj] = tmp
                    end
                elseif i == m
                    if j != 1 && j != n
                        tmp = ones(k,k)
                        tmp[k, :] .= 0
                        Mexp[ii, jj] = tmp
                    elseif j == 1
                        tmp = ones(k,k)
                        tmp[:, 1] .= 0
                        tmp[k, :] .= 0
                        Mexp[ii, jj] = tmp
                    elseif j == n
                        tmp = ones(k,k)
                        tmp[k, :] .= 0
                        tmp[:, k] .= 0
                        Mexp[ii, jj] = tmp
                    end
                else
                    if j == 1
                        tmp = ones(k,k)
                        tmp[:, 1] .= 0
                        Mexp[ii, jj] = tmp
                    elseif j == n
                        tmp = ones(k,k)
                        tmp[:, k] .= 0
                        Mexp[ii, jj] = tmp
                    end
                end
            else
                if M[i,j] == 1 || M[i,j] == -1
                    tmp = M[i,j] .* ones(k,k)
                    Mexp[ii, jj] = tmp
                elseif M[i,j] == 2 || M[i,j] == 4
                    tmp = ones(k,k)
                    tmp[:, Int((k+1)/2)] .= M[i,j]
                    Mexp[ii, jj] = tmp
                elseif M[i,j] == 3 || M[i,j] == 5
                    tmp = ones(k,k)
                    tmp[Int((k+1)/2), :] .= M[i,j]
                    Mexp[ii, jj] = tmp
                elseif M[i,j] == 0
                    if M[i,j+1] == 0 && M[i+1,j] != 0 && M[i,j-1] != 0 &&
                    M[i-1,j] != 0
                        # 1
                        tmp = -1 .* ones(k,k)
                        tmp[[1,k], :] .= 0
                        tmp[:, 1] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] == 0 && M[i,j-1] != 0 &&
                    M[i-1,j] != 0
                        # 2
                        tmp = -1 .* ones(k,k)
                        tmp[1, :] .= 0
                        tmp[:, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] != 0 && M[i,j-1] == 0 &&
                    M[i-1,j] != 0
                        # 3
                        tmp = -1 .* ones(k,k)
                        tmp[[1,k], :] .= 0
                        tmp[:, k] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] != 0 && M[i,j-1] != 0 &&
                    M[i-1,j] == 0
                        # 4
                        tmp = -1 .* ones(k,k)
                        tmp[k, :] .= 0
                        tmp[:, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] == 0 && M[i,j-1] != 0 &&
                    M[i-1,j] != 0
                        # 12
                        tmp = -1 .* ones(k,k)
                        tmp[1, :] .= 0
                        tmp[:, 1] .= 0
                        tmp[k, k] = 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] != 0 && M[i,j-1] == 0 &&
                    M[i-1,j] != 0
                        # 13
                        tmp = -1 .* ones(k,k)
                        tmp[[1,k], :] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] != 0 && M[i,j-1] != 0 &&
                    M[i-1,j] == 0
                        # 14
                        tmp = -1 .* ones(k,k)
                        tmp[k, :] .= 0
                        tmp[:, 1] .= 0
                        tmp[1, k] = 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] == 0 && M[i,j-1] == 0 &&
                    M[i-1,j] != 0
                        # 23
                        tmp = -1 .* ones(k,k)
                        tmp[1, :] .= 0
                        tmp[:, k] .= 0
                        tmp[k, 1] = 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] == 0 && M[i,j-1] != 0 &&
                    M[i-1,j] == 0
                        # 24
                        tmp = -1 .* ones(k,k)
                        tmp[:, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] != 0 && M[i,j-1] == 0 &&
                    M[i-1,j] == 0
                        # 34
                        tmp = -1 .* ones(k,k)
                        tmp[k, :] .= 0
                        tmp[:, k] .= 0
                        tmp[1, 1] = 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] == 0 && M[i,j-1] == 0 &&
                    M[i-1,j] != 0
                        # 123
                        tmp = -1 .* ones(k,k)
                        tmp[1, :] .= 0
                        tmp[k, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] == 0 && M[i,j-1] != 0 &&
                    M[i-1,j] == 0
                        # 124
                        tmp = -1 .* ones(k,k)
                        tmp[:, 1] .= 0
                        tmp[[1,k], k] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] != 0 && M[i,j-1] == 0 &&
                    M[i-1,j] == 0
                        # 134
                        tmp = -1 .* ones(k,k)
                        tmp[k, :] .= 0
                        tmp[1, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] == 0 && M[i,j-1] == 0 &&
                    M[i-1,j] == 0
                        # 234
                        tmp = -1 .* ones(k,k)
                        tmp[:, k] .= 0
                        tmp[[1,k], 1] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] == 0 && M[i+1,j] == 0 && M[i,j-1] == 0 &&
                    M[i-1,j] == 0
                        # 1234
                        tmp = -1 .* ones(k,k)
                        tmp[1, [1,k]] .= 0
                        tmp[k, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    elseif M[i,j+1] != 0 && M[i+1,j] != 0 && M[i,j-1] != 0 &&
                    M[i-1,j] != 0
                        tmp = -1 .* ones(k,k)
                        tmp[[1,k], :] .= 0
                        tmp[:, [1,k]] .= 0
                        Mexp[ii, jj] = tmp
                    else
                        error("Invalid configuration of zeros.")
                    end
                else
                    error("Matrix should only contain 0, 1, 2, 3, 4, or 5.")
                end
            end
        end
    end

    return Mexp

end

