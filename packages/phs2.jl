using NearestNeighbors
using SparseArrays
using LinearAlgebra

###########################################################################

function phi(x, y, phs)
    return sqrt.(x .^ 2 + y .^ 2) .^ phs
end

function phi_x(x, y, phs)
    return phs .* x .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2)
end

function phi_y(x, y, phs)
    return phs .* y .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2)
end

function phiHV(x, y, phs, K)
    z = range(phs - 2*(K-1), step = 2, stop = phs);
    return prod(z .^ 2) .* sqrt.(x .^ 2 + y .^ 2) .^ (phs - 2*K);
end

###########################################################################
"""
    getWeights(xy, XY_o, mn, phs, pol, K)

# Description
    Uses polyharmonic splines and polynomials to get weights for
    approximating the derivative of a function based on function
    values.

# Input
    xy   : single spot where you want to get the derivative (2x1)
    XY_o : locations where function values are known (2 columns)
    mn   : derivative to approximate (1x2)
    phs  : exponent in polyharmonic spline RBF (1,3,5,...)
    pol  : highest degree polynomial in the basis (0,1,2,...)
    K    : Hyperviscosity exponent

# Output
    w : array of approximation weights
"""
function getWeights(xy, XY_o, mn, phs, pol, K)

    # Set ell equal to the number of points in XY
    ell = size(XY_o, 1);
    
    # Initialize the shifted nodes:
    XY = zeros(size(XY_o));

    # Shift so xy is at (0,0)
    for i in 1 : ell
        XY[i,:] = XY_o[i,:] .- xy;
    end

    # # normalize by the maximum radius
    # r = sqrt.( XY[:,1] .^ 2 + XY[:,2] .^2 );
    # r = maximum( r );
    # XY = XY ./ r;

    # Calculate the total number of polynomial basis functions
    nPol = Int((pol + 1) * (pol + 2) / 2);

    # Initialize matrices
    P = zeros(ell, nPol);
    A = zeros(ell + nPol, ell + nPol);
    b = zeros(ell + nPol);

    # Fill in  the polynomial matrix
    if (pol > -1)
        P[:,1] .= 1;
    end
    if (pol > 0)
        P[:,2]  = XY[:,1]                     ;
        P[:,3]  = XY[:,2]                     ;
    end
    if (pol > 1)
        P[:,4]  = XY[:,1] .^ 2                ;
        P[:,5]  = XY[:,1]      .* XY[:,2]     ;
        P[:,6]  =                 XY[:,2] .^ 2;
    end
    if (pol > 2)
        P[:,7]  = XY[:,1] .^ 3                ;
        P[:,8]  = XY[:,1] .^ 2 .* XY[:,2]     ;
        P[:,9]  = XY[:,1]      .* XY[:,2] .^ 2;
        P[:,10] =                 XY[:,2] .^ 3;
    end
    if (pol > 3)
        P[:,11] = XY[:,1] .^ 4                ;
        P[:,12] = XY[:,1] .^ 3 .* XY[:,2]     ;
        P[:,13] = XY[:,1] .^ 2 .* XY[:,2] .^ 2;
        P[:,14] = XY[:,1]      .* XY[:,2] .^ 3;
        P[:,15] =                 XY[:,2] .^ 4;
    end
    if (pol > 4)
        P[:,16] = XY[:,1] .^ 5                ;
        P[:,17] = XY[:,1] .^ 4 .* XY[:,2]     ;
        P[:,18] = XY[:,1] .^ 3 .* XY[:,2] .^ 2;
        P[:,19] = XY[:,1] .^ 2 .* XY[:,2] .^ 3;
        P[:,20] = XY[:,1]      .* XY[:,2] .^ 4;
        P[:,21] =                 XY[:,2] .^ 5;
    end

    # Fill in the polyharmonic spline matrix:
    for i in 1:ell
        for j in 1:ell
            A[i,j] = phi(XY[i,1] - XY[j,1], XY[i,2] - XY[j,2], phs);
        end
    end
    A[1:ell, ell+1:end] = P;
    A[ell+1:end, 1:ell] = transpose(P);

    println()
    println(cond(A))

    # First ell elements of the vector b contain the derivative of each RBF
    # basis function evaluated at 0, and there might be nonzero elements
    # after that which account for the derivatives of the polynomials
    if (ell >= nPol) & (phs >= 2*maximum(mn)+1) & (mod(phs,2) == 1)
        if (mn == [0 0])
            b[1:ell] = phi(-XY[:,1], -XY[:,2], phs);
            if pol >= 0
                b[ell+1] = 1.;
            end
        elseif (mn == [1 0])
            b[1:ell] = phi_x(-XY[:,1], -XY[:,2], phs);
            if pol >=  1
                b[ell+2] = 1.;
            end
        elseif (mn == [0 1])
            b[1:ell] = phi_y(-XY[:,1], -XY[:,2], phs);
            if pol >= 1
                b[ell+3] = 1.;
            end
        elseif (mn == [-1 -1])
            if (K < 1)
                error("The exponent for the Laplacian should be at least 1.")
            end
            if (phs >= 2*K+1)
                b[1:ell] = phiHV(-XY[:,1], -XY[:,2], phs, K);
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
    w = A \ b;

    #remove those weights which will be multiplied by 0
    w = w[1 : ell];

    #return the array of weights
    return w

end

###########################################################################

function test_getWeights(mn, phs, pol, K)

    alp = 1.;

    x = alp .* [-5 -3 -1 1 3 5];
    y = alp .* [-5 -3 -1 1 3 5];
    xy = alp .* [-0.22987; 0.159];

    # x = alp .* [-2. -1. 0. 1. 2.];
    # y = alp .* [-2. -1. 0. 1. 2.];
    # xy = [0; 0];

    X = repeat(x, length(y));
    Y = transpose(repeat(y, length(x)));
    XY = hcat(X[:], Y[:]);

    w = getWeights(xy, XY, mn, phs, pol, K);

    # f(x, y) = 1. .* ones(size(x));
    # f_x(x, y) = 0. .* ones(size(x));
    # f_y(x, y) = 0. * ones(size(x));

    f(x, y) = x .* y;
    f_x(x, y) = y;
    f_y(x, y) = x;

    # f(x, y) = x .^ 2 + y .^ 2;
    # f_x(x, y) = 2. .* x;
    # f_y(x, y) = 2. .* y;

    println()
    println(w' * f(XY[:,1], XY[:,2]) - f(xy[1], xy[2]))
    println()

    w = reshape(w, length(y), length(x))

end

###########################################################################
"""
    getDM(; z = -.9:.1:.9, x = -1:.1:1, m = 1,
        phs = 5, pol = 3, stc = 7)

# Description
    Uses polyharmonic splines and polynomials to get a sparse
    differentiation matrix for approximating the derivative of
    a function.

# Input
    z   : locations where you want to approximate the derivative
    x   : locations where function values are known
    m   : derivative you want to approximate (0,1,2,...)
    phs : exponent in the polyharmonic spline RBF (1,3,5,...)
    pol : highest degree polynomial in the basis (0,1,2,...)
    stc : stencil-size

# Output
    W : sparse differentiation matrix
"""
function getDM(; z = -.9:.1:.9, x = -1:.1:1, m = 1,
    phs = 5, pol = 3, stc = 7)
    
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
"""
    getPeriodicDM(; z = (0:pi/10:2*pi)[1:end-1],
        x = (0:pi/10:2*pi)[1:end-1], m = 1,
        phs = 5, pol = 3, stc = 7, period = 2*pi)

# Description
    Uses polyharmonic splines and polynomials to get a sparse
    differentiation matrix for approximating the derivative of
    a periodic function.

# Input
    z      : locations where you want to approximate the derivative
    x      : locations where function values are known
    m      : derivative you want to approximate (0,1,2,...)
    phs    : exponent in the polyharmonic spline RBF (1,3,5,...)
    pol    : highest degree polynomial in the basis (0,1,2,...)
    stc    : stencil-size
    period : period of the data

# Output
    W : sparse differentiation matrix
"""
function getPeriodicDM(; z = (0:pi/10:2*pi)[1:end-1],
    x = (0:pi/10:2*pi)[1:end-1], m = 1,
    phs = 5, pol = 3, stc = 7, period = 2*pi)

    pad = Int64(round((stc-1)/2))

    x = [x[end-pad+1:end].-period; x; x[1:pad].+period]

    W = getDM(z=z, x=x, m=m, phs=phs, pol=pol, stc=stc)

    W[:,pad+1:2*pad] = W[:,pad+1:2*pad] + W[:,end-pad+1:end]

    W[:,end-2*pad+1:end-pad] = W[:,end-2*pad+1:end-pad] + W[:,1:pad]

    W = W[:,pad+1:end-pad]

    return W

end

###########################################################################

function speedtest_getWeights()

    function myTimer(N, Z)
        @time begin
            for i in 1:N
                w = getWeights(z = Z[i]);
            end
        end
    end

    for n in 0:16
        N = 2^n;
        Z = -3 .+ 6*rand(N);
        myTimer(N, Z);
    end

end

###########################################################################
