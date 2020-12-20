
function freeSlipNoFlux!(U, Tx, Ty, Nx, Ny, indA, indB, indC)

    # Start by just extrapolating the density and energy to the ghost nodes:
    U[indC,[1,4]] = 2 .* U[indB,[1,4]] .- U[indA,[1,4]]

    # Get the tangent velocity on the two non-ghost levels:
    uthA = Tx .* U[indA,2] .+ Ty .* U[indA,3]
    uthB = Tx .* U[indB,2] .+ Ty .* U[indB,3]

    # Linear extrapolation of the tangent velocity to the ghost nodes:
    uthC = 2 .* uthB .- uthA
    
    # Get the normal velocity on the outer-most non-ghost level:
    urB = Nx .* U[indB,2] .+ Ny .* U[indB,3]

    # Set ghost nodes values to enforce no-flux condition:
    urC = -urB
    
    # Set the values of u and v on the ghost nodes:
    U[indC,2] = -Ty .* urC .+ Ny .* uthC
    U[indC,3] =  Tx .* urC .- Nx .* uthC

    return U

end

