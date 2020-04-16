#####################################################################

include("../phs2.jl")

#####################################################################
"""
    getAllDMs(c, nodes, phs, pol, stc, K, a)

# Description
    You put in some parameters and this function gives back all of
    the differentiation matrices (DMs) you will need to solve the
    2D acoustic wave equation on the unit disk.

# Input
    c    : wave speed
    nodes: all of the nodes
    phs  : polyharmonic spline exponent
    pol  : highest degree of polynomial to include in the basis
    stc  : stencil size
    K    : hyperviscosity exponent
    a    : hyperviscosity coefficient

# Output
    Wx  : sparse differentiation matrix for d/dx
    Wy  : sparse differentiation matrix for d/dy
    Whv : sparse differentiation matrix for hyperviscosity operator
"""
function getAllDMs(c, nodes, phs, pol, stc, K, a)

    Wx = getDM(nodes, nodes, [1 0], phs, pol, stc, K)

    Wy = getDM(nodes, nodes, [0 1], phs, pol, stc, K)

    Whv = getDM(nodes, nodes, [-1 -1], phs, pol, stc, K)

    return c*Wx, c*Wy, a*Whv

end

#####################################################################

function initialCondition!(x, y, U)

    # The initial center of the Gaussian bell function
    x0 = .1
    y0 = -.2

    # Initial condition for rho
    U[:,1] = exp.(-10*((x .- x0) .^ 2 .+ (y .- y0) .^ 2))

    # Note that the other variables, u and v, are initiated as zero.

    return U

end

#####################################################################

function ODEfunction!(t, U, dUdt, cWx, cWy, aWhv, bb)
    
    dUdt[:,1] = cWx * U[:,2] .+ cWy * U[:,3] .+ aWhv * U[:,1]

    dUdt[bb,1] .= 0.

    dUdt[:,2] = cWx * U[:,1] .+ aWhv * U[:,2]

    dUdt[:,3] = cWy * U[:,1] .+ aWhv * U[:,3]

    return dUdt

end

#####################################################################








