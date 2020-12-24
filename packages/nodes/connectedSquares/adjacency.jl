
# Takes as input the matrix of zeros and ones defining the room, and returns
# as output a more detailed adjacency matrix which contains information about
# how the different node locations will be treated in terms of refinement and
# boundary condition enforcement.

function adjacency(z1)

    ###########################################################################

    # The actual nodes

    B = zeros(size(z1))

    for i in 2 : size(z1)[1] - 1
        for j in 2 : size(z1)[2] - 1
            if z1[i,j] == 1
                tmp = hcat(z1[i,j+1], z1[i+1,j], z1[i,j-1], z1[i-1,j])
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
                elseif tmp == [1 0 1 1]
                    B[i,j] = 134
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

    # The ghost nodes

    C = zeros(size(z1))

    # First and last columns
    for i in 1 : size(z1)[1]
        if z1[i,2] == 1
            C[i,1] = -1
        else
            C[i,1] = NaN
        end
        if z1[i,end-1] == 1
            C[i,end] = -3
        else
            C[i,end] = NaN
        end
    end

    # First and last rows
    for j in 1 : size(z1)[2]
        if z1[2,j] == 1
            C[1,j] = -2
        else
            C[1,j] = NaN
        end
        if z1[end-1,j] == 1
            C[end,j] = -4
        else
            C[end,j] = NaN
        end
    end

    # Everywhere else
    for i in 2 : size(z1)[1] - 1
        for j in 2 : size(z1)[2] - 1
            if z1[i,j] == 0
                tmp = hcat(z1[i,j+1], z1[i+1,j], z1[i,j-1], z1[i-1,j])
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
                elseif tmp == [1 0 1 1]
                    C[i,j] = -134
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

    # Combine to get the final adjacency matrix A

    A = zeros(size(z1))

    for i in 1 : size(z1)[1]
        for j in 1 : size(z1)[2]
            if z1[i,j] == 1
                A[i,j] = B[i,j]
            else
                A[i,j] = C[i,j]
            end
        end
    end

    ###########################################################################

    return A

end

