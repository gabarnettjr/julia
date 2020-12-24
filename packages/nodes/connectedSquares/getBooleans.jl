
function getBooleans(A)

    m, n = size(A)

    bool_all = BitArray(zeros(m, n))
    bool_noGhost = BitArray(zeros(m, n))

    nanMats = zeros(m, n, 4)
    count = 0

    for i in 1 : m
        for j in 1 : n
            if !isnan(A[i,j])
                count = count + 1
                bool_all[i,j] = true
                if A[i,j] >= 0
                    bool_noGhost[i,j] = true
                end
            else
                nanMats[i,j,:] .= NaN
            end
        end
    end

    zeroVecs = zeros(Float64, count, 4)

    return bool_all, bool_noGhost, nanMats, zeroVecs

end

