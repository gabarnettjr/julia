
function getIndices(A)

    i_n1 = [];  j_n1 = []
    i_n2 = [];  j_n2 = []
    i_n3 = [];  j_n3 = []
    i_n4 = [];  j_n4 = []

    i_n12 = [];  j_n12 = []
    i_n13 = [];  j_n13 = []
    i_n14 = [];  j_n14 = []
    i_n23 = [];  j_n23 = []
    i_n24 = [];  j_n24 = []
    i_n34 = [];  j_n34 = []

    i_n123 = [];  j_n123 = []
    i_n124 = [];  j_n124 = []
    i_n134 = [];  j_n134 = []
    i_n234 = [];  j_n234 = []

    i_n1234 = [];  j_n1234 = []

    m, n = size(A)

    for i in 1 : m
        for j in 1 : n
            if A[i,j] == -1
                i_n1 = vcat(i_n1, i)
                j_n1 = vcat(j_n1, j)
            elseif A[i,j] == -2
                i_n2 = vcat(i_n2, i)
                j_n2 = vcat(j_n2, j)
            elseif A[i,j] == -3
                i_n3 = vcat(i_n3, i)
                j_n3 = vcat(j_n3, j)
            elseif A[i,j] == -4
                i_n4 = vcat(i_n4, i)
                j_n4 = vcat(j_n4, j)
            elseif A[i,j] == -12
                i_n12 = vcat(i_n12, i)
                j_n12 = vcat(j_n12, j)
            elseif A[i,j] == -13
                i_n13 = vcat(i_n13, i)
                j_n13 = vcat(j_n13, j)
            elseif A[i,j] == -14
                i_n14 = vcat(i_n14, i)
                j_n14 = vcat(j_n14, j)
            elseif A[i,j] == -23
                i_n23 = vcat(i_n23, i)
                j_n23 = vcat(j_n23, j)
            elseif A[i,j] == -24
                i_n24 = vcat(i_n24, i)
                j_n24 = vcat(j_n24, j)
            elseif A[i,j] == -34
                i_n34 = vcat(i_n34, i)
                j_n34 = vcat(j_n34, j)
            elseif A[i,j] == -123
                i_n123 = vcat(i_n123, i)
                j_n123 = vcat(j_n123, j)
            elseif A[i,j] == -124
                i_n124 = vcat(i_n124, i)
                j_n124 = vcat(j_n124, j)
            elseif A[i,j] == -134
                i_n134 = vcat(i_n134, i)
                j_n134 = vcat(j_n134, j)
            elseif A[i,j] == -234
                i_n234 = vcat(i_n234, i)
                j_n234 = vcat(j_n234, j)
            elseif A[i,j] == -1234
                i_n1234 = vcat(i_n1234, i)
                j_n1234 = vcat(j_n1234, j)
            end
        end
    end

    return i_n1, j_n1, i_n2, j_n2, i_n3, j_n3, i_n4, j_n4,
           i_n12, j_n12, i_n13, j_n13, i_n14, j_n14, i_n23, j_n23,
           i_n24, j_n24, i_n34, j_n34,
           i_n123, j_n123, i_n124, j_n124, i_n134, j_n134, i_n234, j_n234,
           i_n1234, j_n1234

end

