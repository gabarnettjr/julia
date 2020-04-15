using NearestNeighbors
using SparseArrays
using LinearAlgebra

###########################################################################
#=
The polyharmonic spline radial basis function in 2D, its first derivatives,
and the result of applying the hyperviscosity operator on it K times.
=#
phi(x, y, phs) = sqrt.(x .^ 2 + y .^ 2) .^ phs

phi_x(x, y, phs) = phs .* x .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2)

phi_y(x, y, phs) = phs .* y .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2)

function phiHV(x, y, phs, K)
    z = range(phs - 2*(K-1), step = 2, stop = phs)
    return prod(z .^ 2) .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2*K)
end

###########################################################################
"""
    getWeights(z, x, m, phs, pol, K)

# Description
    Uses polyharmonic splines and polynomials to get weights for
    approximating the derivative of a function based on function
    values.

# Input
    z   : single spot where you want to get the derivative (2x1)
    x   : locations where function values are known (2 rows)
    m   : derivative to approximate (1x2)
    phs : odd exponent in polyharmonic spline RBF (1,3,5,...)
    pol : highest degree polynomial in the basis (0,1,2,3,4,5)
    K   : hyperviscosity exponent (1,2,3,...)

# Output
    w : array of approximation weights
"""
function getWeights(z, x, m, phs, pol, K)

    # Set ell equal to the number of points in x
    ell = size(x, 2)
    
    # Initialize the shifted nodes:
    XY = zeros(size(x))

    # Shift so z is at (0,0)
    for i in 1 : ell
        XY[:,i] = x[:,i] .- z
    end

    # Calculate the total number of 2D polynomial basis functions
    nPol = Int((pol + 1) * (pol + 2) / 2)

    # Initialize matrices
    P = zeros(ell, nPol)
    A = zeros(ell + nPol, ell + nPol)
    b = zeros(ell + nPol)

    # Fill in  the polynomial matrix
    if (pol > -1)
        P[:,1] .= 1
    end
    if (pol > 0)
        P[:,2]  = XY[1,:]
        P[:,3]  = XY[2,:]
    end
    if (pol > 1)
        P[:,4]  = XY[1,:] .^ 2
        P[:,5]  = XY[1,:]      .* XY[2,:]
        P[:,6]  =                 XY[2,:] .^ 2
    end
    if (pol > 2)
        P[:,7]  = XY[1,:] .^ 3
        P[:,8]  = XY[1,:] .^ 2 .* XY[2,:]
        P[:,9]  = XY[1,:]      .* XY[2,:] .^ 2
        P[:,10] =                 XY[2,:] .^ 3
    end
    if (pol > 3)
        P[:,11] = XY[1,:] .^ 4
        P[:,12] = XY[1,:] .^ 3 .* XY[2,:]
        P[:,13] = XY[1,:] .^ 2 .* XY[2,:] .^ 2
        P[:,14] = XY[1,:]      .* XY[2,:] .^ 3
        P[:,15] =                 XY[2,:] .^ 4
    end
    if (pol > 4)
        P[:,16] = XY[1,:] .^ 5
        P[:,17] = XY[1,:] .^ 4 .* XY[2,:]
        P[:,18] = XY[1,:] .^ 3 .* XY[2,:] .^ 2
        P[:,19] = XY[1,:] .^ 2 .* XY[2,:] .^ 3
        P[:,20] = XY[1,:]      .* XY[2,:] .^ 4
        P[:,21] =                 XY[2,:] .^ 5
    end

    # Fill in the polyharmonic spline portion of the matrix A:
    for i in 1:ell
        for j in 1:ell
            A[i,j] = phi(XY[1,i] - XY[1,j], XY[2,i] - XY[2,j], phs)
        end
    end
    A[1:ell, ell+1:end] = P
    A[ell+1:end, 1:ell] = transpose(P)

    # println()
    # println(cond(A))

    # First ell elements of the vector b contain the derivative of each RBF
    # basis function evaluated at 0, and there might be nonzero elements
    # after that which account for the derivatives of the polynomials
    if (ell >= nPol) & (phs >= 2*maximum(m)+1) & (mod(phs,2) == 1)
        if m == [0 0]
            b[1:ell] = phi(-XY[1,:], -XY[2,:], phs)
            if pol > -1
                b[ell+1] = 1.
            end
        elseif m == [1 0]
            b[1:ell] = phi_x(-XY[1,:], -XY[2,:], phs)
            if pol >  0
                b[ell+2] = 1.
            end
        elseif m == [0 1]
            b[1:ell] = phi_y(-XY[1,:], -XY[2,:], phs)
            if pol > 0
                b[ell+3] = 1.
            end
        elseif m == [-1 -1]
            if K < 1
                error("The exponent for the Laplacian should be at least 1.")
            end
            if phs >= 2*K+1
                b[1:ell] = phiHV(-XY[1,:], -XY[2,:], phs, K)
            else
                error("Need phs to be larger to approximate this derivative.")
            end
        else
            error("interpolation, first derivatives, or HV only please.")
        end
    else
        error("Bad parameter values.  Please make sure that phs is an ",
              "odd number, and that length(XY)>=nPol and phs>=2*max(mn)+1.",
              "Also, it could be that the derivative you are requesting",
              "is not yet implemented.")
    end

    #solve linear system for the weights
    w = A \ b

    #remove weights that will be multiplied by 0, and return the rest
    return w[1:ell]

end

###########################################################################

function test_getWeights(m, phs, pol, K)

    alp = 1.

    # x1 = alp .* [-5 -3 -1 1 3 5]
    # x2 = alp .* [-5 -3 -1 1 3 5]
    # z = alp .* [-0.22987; 0.159]

    x1 = alp .* [-2. -1. 0. 1. 2.]
    x2 = alp .* [-2. -1. 0. 1. 2.]
    z = alp .* [0.234; -0.764]

    X = repeat(x1, length(x2), 1)
    Y = repeat(x2', 1, length(x1))
    x = vcat(X[:]', Y[:]')

    w = getWeights(z, x, m, phs, pol, K)

    # f(x, y) = 1. .* ones(size(x))
    # f_x(x, y) = 0. .* ones(size(x))
    # f_y(x, y) = 0. * ones(size(x))

    f(x, y) = x .* y
    f_x(x, y) = y
    f_y(x, y) = x

    # f(x, y) = x .^ 2 + y .^ 2
    # f_x(x, y) = 2. .* x
    # f_y(x, y) = 2. .* y

    println()
    println(dot(w, f(x[1,:], x[2,:])) - f_y(z[1], z[2]))
    println()

    w = reshape(w, length(x2), length(x1))

end

###########################################################################

function speedtest_getWeights()

    x1 = [-5 -3 -1 1 3 5]
    x2 = [-5 -3 -1 1 3 5]

    X = repeat(x1, length(x2), 1)
    Y = repeat(x2', 1, length(x1))
    x = vcat(X[:]', Y[:]')

    m = [0 0]
    phs = 5
    pol = 4
    K = 0

    function myTimer(N, z)
        @time begin
            for i in 1:N
                w = getWeights(z[:,i], x, m, phs, pol, K)
            end
        end
    end

    for n in 0:10
        N = 2^n
        z = -.5 .+ rand(2, N)
        myTimer(N, z)
    end

end

###########################################################################
"""
    getDM(z, x, m, phs, pol, stc, K)

# Description
    Uses polyharmonic splines and polynomials to get a sparse
    differentiation matrix for approximating the derivative of
    a function.

# Input
    z   : locations where you want to approximate the derivative (2 rows)
    x   : locations where function values are known (2 rows)
    m   : derivative you want to approximate (1x2)
    phs : exponent in the polyharmonic spline RBF (1,3,5,...)
    pol : highest degree polynomial in the basis (0,1,2,...)
    stc : stencil-size
    K   : hyperviscosity exponent

# Output
    W : sparse differentiation matrix
"""
function getDM(z, x, m, phs, pol, stc, K)
    
    Lz = size(z, 2)
    Lx = size(x, 2)

    ii = repeat(transpose(1:Lz), stc, 1)
    jj = zeros(Int, stc, Lz)

    w = zeros(stc, Lz)

    tree = KDTree(x)
    idx = knn(tree, z, stc)[1]

    for i in 1:Lz
        jj[:,i] = idx[i]
        w[:,i] = getWeights(z[:,i], x[:,idx[i]], m, phs, pol, K)
    end

    return sparse(ii[:], jj[:], w[:], Lz, Lx)

end

###########################################################################

function test_getDM()

    M = 26
    N = 26

    xx = range(-1, stop=1, length=N)
    yy = range(-1, stop=1, length=M)

    X = repeat(xx', length(yy), 1)
    Y = repeat(yy, 1, length(xx))

    z = vcat(X[:]', Y[:]')
    x = copy(z)
    m = [1 0]
    phs = 5
    pol = 2
    stc = 25
    K = 0

    W = getDM(z, x, m, phs, pol, stc, K)

    f(x,y) = cos.(pi*x) .* sin.(pi*y)
    f_x(x,y) = -pi .* sin.(pi*x) .* sin.(pi*y)
    f_y(x,y) = cos.(pi*x) .* pi .* cos.(pi*y)

    approx = W * f(x[1,:], x[2,:])
    approx = reshape(approx, M, N)

    println(maximum(abs(approx - f_x(X,Y))))

end

###########################################################################
# """
#     getPeriodicDM(; z = (0:pi/10:2*pi)[1:end-1],
#         x = (0:pi/10:2*pi)[1:end-1], m = 1,
#         phs = 5, pol = 3, stc = 7, period = 2*pi)
# 
# # Description
#     Uses polyharmonic splines and polynomials to get a sparse
#     differentiation matrix for approximating the derivative of
#     a periodic function.
# 
# # Input
#     z      : locations where you want to approximate the derivative
#     x      : locations where function values are known
#     m      : derivative you want to approximate (0,1,2,...)
#     phs    : exponent in the polyharmonic spline RBF (1,3,5,...)
#     pol    : highest degree polynomial in the basis (0,1,2,...)
#     stc    : stencil-size
#     period : period of the data
# 
# # Output
#     W : sparse differentiation matrix
# """
# function getPeriodicDM(; z = (0:pi/10:2*pi)[1:end-1],
#     x = (0:pi/10:2*pi)[1:end-1], m = 1,
#     phs = 5, pol = 3, stc = 7, period = 2*pi)
# 
#     pad = Int64(round((stc-1)/2))
# 
#     x = [x[end-pad+1:end].-period; x; x[1:pad].+period]
# 
#     W = getDM(z=z, x=x, m=m, phs=phs, pol=pol, stc=stc)
# 
#     W[:,pad+1:2*pad] = W[:,pad+1:2*pad] + W[:,end-pad+1:end]
# 
#     W[:,end-2*pad+1:end-pad] = W[:,end-2*pad+1:end-pad] + W[:,1:pad]
# 
#     W = W[:,pad+1:end-pad]
# 
#     return W
# 
# end

###########################################################################

