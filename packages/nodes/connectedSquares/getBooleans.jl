
function getBooleans(A)

    m, n = size(A)

    bool_all = BitArray(zeros(m, n))
    bool_noGhost = BitArray(zeros(m, n))

    for i in 1 : m
        for j in 1 : n
            if !isnan(A[i,j])
                bool_all[i,j] = true
                if A[i,j] >= 0
                    bool_noGhost[i,j] = true
                end
            end
        end
    end

    return bool_all, bool_noGhost

end

