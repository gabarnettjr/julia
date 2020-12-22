
# include("zerosAndOnes.jl")

function adjacency(A)

    ###########################################################################

    # the actual nodes

    B = zeros(size(A))

    for i in 2 : size(A)[1] - 1
        for j in 2 : size(A)[2] - 1
            if A[i,j] == 1
                tmp = hcat(A[i,j+1], A[i+1,j], A[i,j-1], A[i-1,j])
                tmp = -1 .* (tmp .- 1)
                if tmp == [1 0 0 0]
                    B[i,j] = 1
                elseif tmp == [0 1 0 0]
                    B[i,j] = 2
                elseif tmp == [0 0 1 0]
                    B[i,j] = 3
                elseif tmp == [0 0 0 1]
                    B[i,j] = 4
                elseif tmp == [1 1 0 0]
                    B[i,j] = 12
                elseif tmp == [1 0 1 0]
                    B[i,j] = 13
                elseif tmp == [1 0 0 1]
                    B[i,j] = 14
                elseif tmp == [0 1 1 0]
                    B[i,j] = 23
                elseif tmp == [0 1 0 1]
                    B[i,j] = 24
                elseif tmp == [0 0 1 1]
                    B[i,j] = 34
                elseif tmp == [1 1 1 0]
                    B[i,j] = 123
                elseif tmp == [1 1 0 1]
                    B[i,j] = 124
                elseif tmp == [0 1 1 1]
                    B[i,j] = 234
                elseif tmp == [1 1 1 1]
                    B[i,j] = 1234
                else
                    B[i,j] = 0
                end
            end
        end
    end

    ###########################################################################

    # the ghost nodes

    C = zeros(size(A))

    # first and last columns
    for i in 1 : size(A)[1]
        if A[i,2] == 1
            C[i,1] = -1
        else
            C[i,1] = NaN
        end
        if A[i,end-1] == 1
            C[i,end] = -3
        else
            C[i,end] = NaN
        end
    end

    # first and last rows
    for j in 1 : size(A)[2]
        if A[2,j] == 1
            C[1,j] = -2
        else
            C[1,j] = NaN
        end
        if A[end-1,j] == 1
            C[end,j] = -4
        else
            C[end,j] = NaN
        end
    end

    # everywhere else
    for i in 2 : size(A)[1] - 1
        for j in 2 : size(A)[2] - 1
            if A[i,j] == 0
                tmp = hcat(A[i,j+1], A[i+1,j], A[i,j-1], A[i-1,j])
                if tmp == [1 0 0 0]
                    C[i,j] = -1
                elseif tmp == [0 1 0 0]
                    C[i,j] = -2
                elseif tmp == [0 0 1 0]
                    C[i,j] = -3
                elseif tmp == [0 0 0 1]
                    C[i,j] = -4
                elseif tmp == [1 1 0 0]
                    C[i,j] = -12
                elseif tmp == [1 0 1 0]
                    C[i,j] = -13
                elseif tmp == [1 0 0 1]
                    C[i,j] = -14
                elseif tmp == [0 1 1 0]
                    C[i,j] = -23
                elseif tmp == [0 1 0 1]
                    C[i,j] = -24
                elseif tmp == [0 0 1 1]
                    C[i,j] = -34
                elseif tmp == [1 1 1 0]
                    C[i,j] = -123
                elseif tmp == [1 1 0 1]
                    C[i,j] = -124
                elseif tmp == [0 1 1 1]
                    C[i,j] = -234
                elseif tmp == [1 1 1 1]
                    C[i,j] = -1234
                else
                    C[i,j] = NaN
                end
            end
        end
    end

    ###########################################################################

    # combine to get the final adjacency matrix D

    D = zeros(size(A))

    for i in 1 : size(A)[1]
        for j in 1 : size(A)[2]
            if A[i,j] == 1
                D[i,j] = B[i,j]
            else
                D[i,j] = C[i,j]
            end
        end
    end

    ###########################################################################

    return D

end

