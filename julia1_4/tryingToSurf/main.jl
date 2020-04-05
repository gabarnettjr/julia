using Plots

gr()

M = 25
N = 50

region = "annulus"

keepGoing = true

function myMeshgrid( x, y )
    xx = repeat( x', length(y) )
    yy = repeat( y', length(x) )
    return xx, yy'
end

if region == "rectangle"

    x = range( -1, 1, length = N+1 )
    y = range( -1, 1, length = M+1 )
    
    xx, yy = myMeshgrid( x, y )

elseif region == "annulus"

    r = range( 1, 2, length = N+1 )
    th = range( 0, 2*pi, length = M+1 )
    th = th[ 1 : end-1 ]
    
    rr, thth = myMeshgrid( r, th )
    
    xx = rr .* cos.(thth)
    yy = rr .* sin.(thth)

else

    println()
    println("Please choose rectangle or annulus for the region string.")
    keepGoing = false

end

if keepGoing

    X = xx[:]
    Y = yy[:]
    
    Z = sqrt.( X .^ 2 .+ Y .^ 2 )
    # Z = ones( size(X) )
    
    surface( X, Y, Z, camera = (0,90), grid = false )

end




