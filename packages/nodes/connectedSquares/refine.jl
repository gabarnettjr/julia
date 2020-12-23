
function refine!(z1, A, refinement)
    
    if refinement == 1
        return z1, A
    end

    m = size(A)[1]
    n = size(A)[2]

    z1_new = zeros(2*m, 2*n)
    A_new = zeros(2*m, 2*n)
    
    for k = 1 : (refinement - 1)

        for i in 1 : m

            for j in 1 : n

                ii = 2*i-1 : 2*i
                jj = 2*j-1 : 2*j

                if z1[i,j] == 1
                    z1_new[ii, jj] = ones(2,2)
                else
                    z1_new[ii, jj] = zeros(2,2)
                end
                
                # The actual squares
                if A[i,j] == 1
                    A_new[ii,jj] = [0 1; 0 1]
                elseif A[i,j] == 2
                    A_new[ii,jj] = [2 2; 0 0]
                elseif A[i,j] == 3
                    A_new[ii,jj] = [3 0; 3 0]
                elseif A[i,j] == 4
                    A_new[ii,jj] = [0 0; 4 4]
                elseif A[i,j] == 12
                    A_new[ii,jj] = [2 12; 0 1]
                elseif A[i,j] == 13
                    A_new[ii,jj] = [3 1; 3 1]
                elseif A[i,j] == 14
                    A_new[ii,jj] = [0 1; 4 14]
                elseif A[i,j] == 23
                    A_new[ii,jj] = [23 2; 3 0]
                elseif A[i,j] == 24
                    A_new[ii,jj] = [2 2; 4 4]
                elseif A[i,j] == 34
                    A_new[ii,jj] = [3 0; 34 4]
                elseif A[i,j] == 123
                    A_new[ii,jj] = [23 12; 3 1]
                elseif A[i,j] == 124
                    A_new[ii,jj] = [2 12; 4 14]
                elseif A[i,j] == 134
                    A_new[ii,jj] = [3 1; 34 4]
                elseif A[i,j] == 234
                    A_new[ii,jj] = [23 2; 34 4]
                elseif A[i,j] == 1234
                    A_new[ii,jj] = [23 12; 34 14]
                end

                # The ghost squares
                if A[i,j] == -1
                    A_new[ii,jj] = [NaN -1; NaN -1]
                elseif A[i,j] == -2
                    A_new[ii,jj] = [-2 -2; NaN NaN]
                elseif A[i,j] == -3
                    A_new[ii,jj] = [-3 NaN; -3 NaN]
                elseif A[i,j] == -4
                    A_new[ii,jj] = [NaN NaN; -4 -4]
                elseif A[i,j] == -12
                    A_new[ii,jj] = [-2 -12; NaN -1]
                elseif A[i,j] == -13
                    A_new[ii,jj] = [-3 -1; -3 -1]
                elseif A[i,j] == -14
                    A_new[ii,jj] = [NaN -1; -4 -14]
                elseif A[i,j] == -23
                    A_new[ii,jj] = [-23 -2; -3 NaN]
                elseif A[i,j] == -24
                    A_new[ii,jj] = [-2 -2; -4 -4]
                elseif A[i,j] == -34
                    A_new[ii,jj] = [-3 NaN; -34 -4]
                elseif A[i,j] == -123
                    A_new[ii,jj] = [-23 -12; -3 -1]
                elseif A[i,j] == -124
                    A_new[ii,jj] = [-2 -12; -4 -14]
                elseif A[i,j] == -134
                    A_new[ii,jj] = [-3 -1; -34 -4]
                elseif A[i,j] == -234
                    A_new[ii,jj] = [-23 -2; -34 -4]
                elseif A[i,j] == -1234
                    A_new[ii,jj] = [-23 -12; -34 -14]
                elseif A[i,j] == NaN
                    A_new[ii,jj] = [NaN NaN; NaN NaN]
                end

            end

        end

    end

    return z1_new, A_new

end

