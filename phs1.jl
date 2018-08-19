using NearestNeighbors
using SparseArrays

###########################################################################

function getWeights(; z=0, x=-3:3, m=1, phs=5, pol=3)
    #=
    Uses polyharmonic splines and polynomials to get weights for
    approximating the derivative of a function based on function
    values.

    INPUT
    z   : location where you want to approximate the derivative
    x   : locations where function values are known
    m   : derivative to approximate (0,1,2,...)
    phs : exponent in polyharmonic spline RBF (1,3,5,...)
    pol : highest degree polynomial in the basis (0,1,2,...)

    OUTPUT
    w : array of approximation weights
    =#

    width = maximum(x) - minimum(x)                      #width of interval
    ell = length(x)                                     #length of vector x

    x = (x .- z) / width      #shift to zero and scale by width of interval

    #Initialize matrices:
    P = zeros(ell, pol+1)
    A = zeros(ell+pol+1, ell+pol+1)
    b = zeros(ell+pol+1)

    #Fill in polynomial matrix:
    for j in 0:pol
        P[:,j+1] = x.^j
    end

    #Fill in full polyharmonic spline plus polynomial matrix:
    for i in 1:ell
        for j in 1:ell
            A[i,j] = abs.(x[i] - x[j]) .^ phs
        end
    end
    A[1:ell, ell+1:end] = P
    A[ell+1:end, 1:ell] = transpose(P)

    #First ell elements of the vector b contain the derivative of each RBF
    #basis function evaluated at 0:
    if (ell >= pol+1) & (phs >= m+1) & (mod(phs,2) == 1)
        if mod(m,2) == 0
            b[1:ell] = prod(phs-(m-1) : phs) .*
                abs.(-x) .^ (phs-m)
        else
            b[1:ell] = prod(phs-(m-1) : phs) .*
                (-x) .^ (phs-m) .* sign.(-x)
        end
    else
        error("Bad parameter values.  Please make sure that phs is an ",
            "odd number, and that length(x)>=pol+1 and phs>=m+1.")
    end

    #Last elements of vector b contain the derivative of each monomial
    #basis function evaluated at 0 (these derivatives are all zero except
    #maybe one of them):
    if pol >= m
        b[ell+m+1] = factorial(m)
    end

    w = A \ b                          #solve linear system for the weights

    w = w[1:ell]        #remove those weights which will be multiplied by 0

    return w ./ width^m                    #return correctly scaled weights

end

###########################################################################

function getDM(; z=-.9:.1:.9, x=-1:.1:1, m=1,
    phs=5, pol=3, stc=7)
    #=
    Uses polyharmonic splines and polynomials to get a sparse
    differentiation matrix for approximating the derivative of
    a function.

    INPUT
    z   : locations where you want to approximate the derivative
    x   : locations where function values are known
    m   : derivative you want to approximate (0,1,2,...)
    phs : exponent in the polyharmonic spline RBF (1,3,5,...)
    pol : highest degree polynomial in the basis (0,1,2,...)
    stc : stencil-size

    OUTPUT
    W : sparse differentiation matrix
    =#
    
    Lx = length(x)
    Lz = length(z)

    tmp = zeros(2, Lx)
    tmp[1,:] = x
    x = tmp

    tmp = zeros(2, Lz)
    tmp[1,:] = z
    z = tmp

    ii = repeat(transpose(1:Lz), stc, 1)
    jj = zeros(Int64, stc, Lz)

    w = zeros(stc, Lz)

    tree = KDTree(x)

    for i in 1:Lz
        jj[:,i] = knn(tree, z[:,i], stc)[1]
        w[:,i] = getWeights(z=z[1,i], x=x[1,jj[:,i]], m=m,
            phs=phs, pol=pol)
    end

    return sparse(ii[:], jj[:], w[:], Lz, Lx)

end

###########################################################################

function getPeriodicDM(; z=(0:pi/10:2*pi)[1:end-1],
    x=(0:pi/10:2*pi)[1:end-1], m=1,
    phs=5, pol=3, stc=7, period=2*pi)
    #=
    Uses polyharmonic splines and polynomials to get a sparse
    differentiation matrix for approximating the derivative of
    a periodic function.

    INPUT
    z      : locations where you want to approximate the derivative
    x      : locations where function values are known
    m      : derivative you want to approximate (0,1,2,...)
    phs    : exponent in the polyharmonic spline RBF (1,3,5,...)
    pol    : highest degree polynomial in the basis (0,1,2,...)
    stc    : stencil-size
    period : period of the data

    OUTPUT
    W : sparse differentiation matrix
    =#

    pad = Int64(round((stc-1)/2))

    x = [x[end-pad+1:end].-period; x; x[1:pad].+period]

    W = getDM(z=z, x=x, m=m, phs=phs, pol=pol, stc=stc)

    W[:,pad+1:2*pad] = W[:,pad+1:2*pad] + W[:,end-pad+1:end]

    W[:,end-2*pad+1:end-pad] = W[:,end-2*pad+1:end-pad] + W[:,1:pad]

    W = W[:,pad+1:end-pad]

    return W

end

###########################################################################
