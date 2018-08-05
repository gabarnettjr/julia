function getWeights( z, x, m, phs, pol )

    width = maximum(x) - minimum(x)                      #width of interval
    ell = length(x)                                     #length of vector x

    x = ( x - z ) / width     #shift to zero and scale by width of interval

    #Initialize matrices if necessary:
    P = zeros( ell, pol+1 )
    A = zeros( ell+pol+1, ell+pol+1 )
    b = zeros( ell+pol+1 )

    #Make polynomial matrix:
    for j in 0:pol
        P[:,j+1] = x.^j
    end

    #Make the full polyharmonic spline plus polynomial matrix:
    for i in 1:ell
        for j in 1:ell
            A[i,j] = abs( x[i] - x[j] ) .^ phs
        end
    end
    A[ 1:ell, ell+1:end ] = P
    A[ ell+1:end, 1:ell ] = P.'

    #First ell elements of the vector b contain the derivative of each RBF
    #basis function evaluated at 0:
    if ( ell >= pol+1 ) & ( phs >= m+1 )
        if mod(m,2) == 0
            b[1:ell] = prod( phs-(m-1) : phs ) .*
                abs(0.-x) .^ (phs-m)
        else
            b[1:ell] = prod( phs-(m-1) : phs ) .*
                (0.-x) .^ (phs-m) .* sign.(collect(0.-x))
        end
    else
        error( "Bad parameter values." )
    end

    #Last elements of vector b contain the derivative of each monomial
    #basis function evaluated at 0 (these derivatives are all zero except
    #maybe one of them):
    if pol >= m
        b[ell+m+1] = factorial(m)
    end

    w = A \ b                          #solve linear system for the weights

    w = w[1:ell]          #remove the weights which will be multiplied by 0

    return w ./ width^m                    #return correctly scaled weights

end

###########################################################################

using NearestNeighbors

function getDM( x, X; m=1, phs=5, pol=3, stc=7 )
    
    ell = length(x)
    L = length(X)

    tmp = zeros( 2, ell )
    tmp[1,:] = x
    x = tmp

    tmp = zeros( 2, L )
    tmp[1,:] = X
    X = tmp

    ii = repmat( (1:L).', stc, 1 )
    jj = zeros( Int64, size(ii) )

    w = zeros( stc, L )

    tree = KDTree( x )

    for i in 1:L
        ( jj[:,i], d ) = knn( tree, X[:,i], stc )
        w[:,i] = getWeights( X[1,i], x[1,jj[:,i]], m, phs, pol )
    end

    return sparse( ii[:], jj[:], w[:], L, ell )

end

###########################################################################

function getPeriodicDM( x, X; m=1, phs=5, pol=3, stc=7, period=2*pi )

    pad = Int64( round( (stc+1)/2 ) )

    x = [ x[end-pad:end]-period, x, x[1:pad]+period ]

    ell = length(x)

    W = getDM( x, X, m, phs, pol, stc )

    return W     #NOT FINISHED YET

end

###########################################################################
