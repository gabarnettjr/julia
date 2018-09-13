using Printf
using PyPlot
using LinearAlgebra
using DelimitedFiles

include("../rk.jl")
include("../phs1.jl")

###########################################################################

testCase = "risingBubble"

verticallyLagrangian = false

#Choose 0, 1, 2, 3, or 4:
refinementLevel = 1

saveArrays       = false
saveContours     = false
contourFromSaved = true
plotNodesAndExit = false

#Choose "u", "w", "T", "rho", "phi", "P", "theta", "pi", or "phi":
whatToPlot = "w"

#Choose either a number of contours, or a range of contours:
contours = 20

###########################################################################

#Get string for saving results:

saveString = string("./results/", testCase, "_")

if verticallyLagrangian
    saveString = string(saveString, "vLag", "_")
else
    saveString = string(saveString, "vEul", "_")
end

saveString = string(saveString, refinementLevel, "/")

###########################################################################k

if ( saveArrays ) & ( ! contourFromSaved )
    if ! isdir(saveString)
        mkpath(saveString)                     #make new directories
    end
end

if saveContours
    tmp = walkdir("results")
    for (roots, dirs, files) in tmp
        for file in files
            if endswith(file, ".png") | endswith(file, ".PNG")
                if Sys.islinux()
                    # print(string("nonhydroXZresults","/",file,"\n"))
                    rm(string("results","/",file))
                elseif Sys.iswindows()
                    # print(string("nonhydroXZresults","\\",file,"\n"))
                    rm(string("results","\\",file))
                else
                    error("Only implemented for linux and windows.")
                end
            end
        end
    end
end

if contourFromSaved
    saveArrays = false
    saveContours = true
end

###########################################################################

#Definitions of constants:

Cp = 1004                        #specific heat of air at constant pressure
Cv = 717                           #specific heat of air at constant volume
Rd = Cp - Cv                                      #gas constant for dry air
g  = 9.81                                           #gravitational constant
Po = 10^5                                    #reference pressure in Pascals

th0 = 300                                  #reference potential temperature
N   = .01                                          #Brunt-Vaisala frequency

###########################################################################

#All of the test cases are defined in terms of hydrostatic background
#states for potential temperature (th) and Exner pressure (pi), and these
#are functions of z only.  In the pressure coordinate case, the initial
#s-levels are defined using the inverse Exner pressure function.

if (testCase == "risingBubble") | (testCase == "densityCurrent")

    function potentialTemperature(z)
        return th0 * ones(size(z))
    end

    function potentialTemperatureDerivative(z)
        return zeros(size(z))
    end

    function exnerPressure(z)
        return 1 .- g/Cp/th0 * z
    end

    function inverseExnerPressure(pEx)
        return (1 .- pEx) * Cp*th0/g
    end

elseif testCase == "inertiaGravityWaves"

    function potentialTemperature(z)
        return th0 * exp(N^2/g * z)
    end

    function potentialTemperatureDerivative(z)
        return th0 * N^2/g * exp(N^2/g * z)
    end

    function exnerPressure(z)
        return 1 + g^2/Cp/th0/N^2 * (exp(-N^2/g*z) .- 1)
    end

    function inverseExnerPressure(pi)
        return -g/N^2 * log(1 + (pi-1) * Cp*th0*N^2/g^2)
    end

else

    error("Invalid test case string.")

end

###########################################################################

#Get some test-specific parameters, such as the size of the domain, the
#horizontal node spacing dx, and the number of vertical levels nLev:

if testCase == "risingBubble"
    left    = -5000
    right   = 5000
    bottom  = 0
    top     = 10000
    dx      = 400 / 2^refinementLevel
    nLev    = 25 * 2^refinementLevel
    dt      = 1/2 / 2^refinementLevel
    tf      = 1200
    saveDel = 100
    function zSurfFunction(x)
        return 1000 * sin.(2*pi/10000 * x)
    end
    function zSurfPrimeFunction(x)
        return 1000 * 2*pi/10000 * cos.(2*pi/10000 * x)
    end
    function thetaPtb(x, z)
        return 2 * exp.(-1e-6*(x.^2 + (z .- 4000).^2))
    end
elseif testCase == "densityCurrent"
    left    = -25600
    right   = 25600
    bottom  = 0
    top     = 6400
    dx      = 400 / 2^refinementLevel
    nLev    = 16 * 2^refinementLevel
    dt      = 1/3 / 2^refinementLevel
    tf      = 900
    saveDel = 50
    function zSurfFunction(x)
        return zeros(size(x))
    end
    function zSurfPrimeFunction(x)
        return zeros(size(x))
    end
    function thetaPtb(x, z)
        return -20 * exp(-7e-7*(x.^2 + (z-3000).^2))
    end
elseif testCase == "inertiaGravityWaves"
    left    = -150000
    right   = 150000
    bottom  = 0
    top     = 10000
    dx      = 1000 / 2^refinementLevel
    nLev    = 10 * 2^refinementLevel
    dt      = 1/1 / 2^refinementLevel
    tf      = 3000
    saveDel = 250
    function zSurfFunction(x)
        return zeros(size(x))
    end
    function zSurfPrimeFunction(x)
        return zeros(size(x))
    end
    function thetaPtb(x, z)
        thetaC = .01
        hC = 10000
        aC = 5000
        xC = -50000
        return thetaC * sin.(pi * z / hC) ./ (1 .+ ((x - xC) / aC).^2)
    end
else
    error("Invalid test case string.")
end

###########################################################################

t = 0

nTimesteps = Int64(round(tf / dt))

nCol = Int64(round((right - left) / dx))

x = range(left+dx/2, stop=right-dx/2, length=nCol)

zSurf = zSurfFunction(x)
zSurfPrime = zSurfPrimeFunction(x)

###########################################################################

function getSvalues()

    ds = (top - bottom) / nLev
    s = range(bottom-ds/2, stop=top+ds/2, length=nLev+2)

    return s, ds

end

###########################################################################

s, ds = getSvalues()

xx = transpose(repeat(x, 1, nLev+2))
ss = repeat(s, 1, nCol)

###########################################################################

#All of the polyharmonic spline radial basis function weights:

phs = 11                                 #lateral PHS exponent (odd number)
pol = 5                         #highest degree polynomial in lateral basis
stc = 11                                              #lateral stencil size
alp = 2^-14 * 300                          #lateral dissipation coefficient
Wa   = getPeriodicDM(z=x, x=x, m=1,
    phs=phs, pol=pol, stc=stc, period=right-left)
Whva = getPeriodicDM(z=x, x=x, m=pol+1,
    phs=phs, pol=pol, stc=stc, period=right-left)
Whva = alp * dx^pol * Whva

Wa = copy(transpose(Wa))
Whva = copy(transpose(Whva))

phs = 5                                              #vertical PHS exponent
pol = 3                        #highest degree polynomial in vertical basis
stc = 7                                              #vertical stencil size
alp = -2^-10 * 300                        #vertical dissipation coefficient
Ws   = getDM(z=s,          x=s, m=1,     phs=phs, pol=pol, stc=stc)
Whvs = getDM(z=s[2:end-1], x=s, m=pol+1, phs=phs, pol=pol, stc=stc)
Whvs = alp * ds^pol * Whvs

wIbot = getWeights(z=(s[1]+s[2])/2.,   x=s[1:stc],        m=0,
    phs=phs, pol=pol)
wEbot = getWeights(z=s[1],             x=s[2:stc+1],      m=0,
    phs=phs, pol=pol)
wDbot = getWeights(z=(s[1]+s[2])/2.,   x=s[1:stc],        m=1,
    phs=phs, pol=pol)
wHbot = getWeights(z=(s[1]+s[2])/2.,   x=s[2:stc+1],      m=0,
    phs=phs, pol=pol)

wItop = getWeights(z=(s[end-1]+s[end])/2., x=s[end:-1:end-stc+1], m=0,
    phs=phs, pol=pol)
wEtop = getWeights(z=s[end],               x=s[end-1:-1:end-stc], m=0,
    phs=phs, pol=pol)
wDtop = getWeights(z=(s[end-1]+s[end])/2., x=s[end:-1:end-stc+1], m=1,
    phs=phs, pol=pol)
wHtop = getWeights(z=s[end],               x=s[end-1:-1:end-stc], m=0,
    phs=phs, pol=pol)

function Da(U)
    return U * Wa
end

function Ds(U)
    return Ws * U
end

function HV(U)
    return U[2:end-1,:] * Whva + Whvs * U
end

###########################################################################

function backgroundStatesAndPerturbations(zz)

    thetaBar = potentialTemperature(zz)
    piBar = exnerPressure(zz)
    piPtb = zeros(Float64, nLev+2, nCol)
    Tbar = piBar .* thetaBar
    Tptb = (piBar + piPtb) .* (thetaBar + thetaPtb(xx,zz)) - Tbar
    Pbar = Po * piBar .^ (Cp/Rd)
    Pptb = Po * (piBar + piPtb) .^ (Cp/Rd) - Pbar
    rhoBar = Pbar / Rd ./ Tbar
    rhoPtb = (Pbar + Pptb) / Rd ./ (Tbar + Tptb) - rhoBar
    phiBar = g * zz

    return thetaBar, piBar, Tbar, Pbar, rhoBar, phiBar,
        Tptb, rhoPtb

end

###########################################################################

function getVerticalLevels()

    zz = zeros(nLev+2, nCol)

    for j in 1:nCol
        dz = (top - zSurf[j]) / (top - bottom) * ds
        zz[:,j] = range(zSurf[j]-dz/2, stop=top+dz/2, length=nLev+2)
    end

    return zz

end

###########################################################################

zz = getVerticalLevels()

if plotNodesAndExit
    figure()
    plot(xx[:], zz[:], marker=".", linestyle="none")
    plot(x, zSurf, linestyle="-", color="red")
    plot(x, top*ones(size(x)), linestyle="-", color="red")
    axis(:image)
    show()
    exit()
end

###########################################################################

#Assignment of hydrostatic background states and initial perturbations:
thetaBar, piBar, Tbar, Pbar, rhoBar, phiBar,
    Tptb, rhoPtb = backgroundStatesAndPerturbations(zz)

#Assignment of initial conditions:
U = zeros(Float64, nLev+2, nCol, 6)
if testCase == "inertiaGravityWaves"
    U[:,:,1] =  20 * ones(Float64, nLev+2, nCol)       #horizontal velocity
end
U[:,:,2] = zeros(Float64, nLev+2, nCol)                  #vertical velocity
U[:,:,3] = Tptb                                                #temperature
U[:,:,4] = rhoPtb                                                  #density
U[:,:,5] = phiBar                                             #geopotential
U[:,:,6] = zeros(Float64, nLev+2, nCol)                           #pressure

###########################################################################

#Unit tangent and unit normal vectors along bottom and top boundaries:

TzBot = zSurfPrime
TxBot = ones(size(TzBot))
tmp = sqrt.(TxBot.^2 + TzBot.^2)
TxBot = TxBot ./ tmp
TzBot = TzBot ./ tmp

NxBot = transpose(repeat(-TzBot, 1, stc-1))
NzBot = transpose(repeat(TxBot, 1, stc-1))

TxBot = transpose(repeat(TxBot, 1, stc))
TzBot = transpose(repeat(TzBot, 1, stc))

NxTop = zeros(Float64, stc-1, nCol)
NzTop = ones(Float64, stc-1, nCol)

TxTop = ones(Float64, stc, nCol)
TzTop = zeros(Float64, stc, nCol)

###########################################################################

if saveContours
    if testCase == "inertiaGravityWaves"
        fig = figure(figsize = (18,3))
    else
        fig = figure(figsize = (18,14))
    end
end

###########################################################################

function contourSomething(U, t)

    if whatToPlot == "u"
        tmp = U[:,:,1]
    elseif whatToPlot == "w"
        tmp = U[:,:,2]
    elseif whatToPlot == "T"
        # tmp = Tbar
        tmp = U[:,:,3]
    elseif whatToPlot == "rho"
        # tmp = rhoBar
        tmp = U[:,:,4]
    elseif whatToPlot == "phi"
        # tmp = phiBar
        tmp = U[:,:,5] - phiBar
    elseif whatToPlot == "P"
        # tmp = Pbar
        tmp = U[:,:,6]
    elseif whatToPlot == "theta"
        # tmp = thetaBar
        tmp = (U[:,:,3]+Tbar) ./ ((U[:,:,6]+Pbar)/Po).^(Rd/Cp) - thetaBar
    elseif whatToPlot == "pi"
        # tmp = piBar
        tmp = ((U[:,:,6]+Pbar) / Po) .^ (Rd/Cp) - piBar
    else
        error("Invalid whatToPlot string.")
    end

    zz = U[:,:,5] ./ g

    clf()
    contourf(xx, zz, tmp, contours)
    if testCase != "inertiaGravityWaves"
        axis(:image)
        colorbar(orientation="vertical")
    else
        colorbar(orientation="horizontal")
    end
    savefig(@sprintf("%04d.png",Int64(round(t))),
        bbox_inches="tight")

end

###########################################################################

V = zeros(Float64, nLev+2, nCol, 4 )

###########################################################################

function verticalRemap(U, z, Z, V)      #used only in vertically Lagrangian
    """
    Interpolate columns of U from z to Z
    nLev is the number of interior levels of U
    V is the output array
    """
    z = repeat(z[:,:,1], 1, 1, size(U,3))
    Z = repeat(Z[:,:,1], 1, 1, size(U,3))
    #quadratic on bottom:
    z0 = z[1,:,:]
    z1 = z[2,:,:]
    z2 = z[3,:,:]
    ZZ = Z[1,:,:]
    V[1,:,:] =
      (ZZ - z1) .* (ZZ - z2) .* U[1,:,:] ./ (z0 - z1) ./ (z0 - z2) +
      (ZZ - z0) .* (ZZ - z2) .* U[2,:,:] ./ (z1 - z0) ./ (z1 - z2) +
      (ZZ - z0) .* (ZZ - z1) .* U[3,:,:] ./ (z2 - z0) ./ (z2 - z1)
    #quadratic on interior:
    z0 = z[1:nLev+0,:,:]
    z1 = z[2:nLev+1,:,:]
    z2 = z[3:nLev+2,:,:]
    ZZ = Z[2:nLev+1,:,:]
    V[2:nLev+1,:,:] =
    (ZZ - z1) .* (ZZ - z2) .* U[1:nLev+0,:,:] ./ (z0 - z1) ./ (z0 - z2) +
    (ZZ - z0) .* (ZZ - z2) .* U[2:nLev+1,:,:] ./ (z1 - z0) ./ (z1 - z2) +
    (ZZ - z0) .* (ZZ - z1) .* U[3:nLev+2,:,:] ./ (z2 - z0) ./ (z2 - z1)
    #quadratic on top:
    z0 = z[nLev+0,:,:]
    z1 = z[nLev+1,:,:]
    z2 = z[nLev+2,:,:]
    ZZ = Z[nLev+2,:,:]
    V[nLev+2,:,:] =
      (ZZ - z1) .* (ZZ - z2) .* U[nLev+0,:,:] ./ (z0 - z1) ./ (z0 - z2) +
      (ZZ - z0) .* (ZZ - z2) .* U[nLev+1,:,:] ./ (z1 - z0) ./ (z1 - z2) +
      (ZZ - z0) .* (ZZ - z1) .* U[nLev+2,:,:] ./ (z2 - z0) ./ (z2 - z1)
    return V
end

###########################################################################

function fastBackgroundStates(phi)

    z = phi / g

    thetaBar = potentialTemperature(z)
    piBar = exnerPressure(z)
    dthetaBarDz = potentialTemperatureDerivative(z)

    Tbar = piBar .* thetaBar
    Pbar = Po * piBar .^ (Cp/Rd)
    rhoBar = Pbar / Rd ./ Tbar

    dpiBarDz = -g ./ Cp ./ thetaBar                  #hydrostatic condition
    dTbarDz = piBar .* dthetaBarDz + thetaBar .* dpiBarDz
    dPbarDz = Po * Cp/Rd * piBar.^(Cp/Rd-1) .* dpiBarDz
    drhoBarDz = (dPbarDz - Rd*rhoBar.*dTbarDz) ./ (Rd * Tbar)

    return Pbar, rhoBar, Tbar, drhoBarDz, dTbarDz

end

###########################################################################

function setGhostNodes(U)

    #Enforce phi=g*z on bottom boundary:
    U[1,:,5] = (g*zSurf' - wIbot[2:stc]' * U[2:stc,:,5]) / wIbot[1]

    #Enforce phi=g*z on top boundary:
    U[end,:,5] = (g*top .-
        wItop[2:stc]' * U[end-1:-1:end-stc+1,:,5]) / wItop[1]

    #Get background states on possibly changing geopotential levels:
    Pbar, rhoBar, Tbar, drhoBarDz, dTbarDz =
        fastBackgroundStates(U[:,:,5])

    #extrapolate tangent velocity uT to bottom ghost nodes:
    uT = U[2:stc+1,:,1] .* TxBot + U[2:stc+1,:,2] .* TzBot
    uT = wEbot' * uT

    #get normal velocity uN on bottom ghost nodes:
    uN = U[2:stc,:,1] .* NxBot + U[2:stc,:,2] .* NzBot
    uN = -wIbot[2:stc]' * uN / wIbot[1]

    #use uT and uN to get (u,w) on bottom ghost nodes:
    U[1,:,1] = uT' .* TxBot[1,:] + uN' .* NxBot[1,:]
    U[1,:,2] = uT' .* TzBot[1,:] + uN' .* NzBot[1,:]

    #get (u,w) on top ghost nodes (easier because it's flat):
    U[end,:,1] = wEtop' * U[end-1:-1:end-stc,:,1]
    U[end,:,2] = -wItop[2:stc]' * U[end-1:-1:end-stc+1,:,2] / wItop[1]

    #get pressure on interior nodes using the equation of state:
    U[2:end-1,:,6] =
        ((rhoBar+U[:,:,4]) * Rd .* (Tbar+U[:,:,3]) - Pbar)[2:end-1,:,1]

    #set pressure on bottom ghost nodes:
    dPda = (wHbot' * U[2:stc+1,:,6]) * Wa
    rho = wHbot' * U[2:stc+1,:,4]
    dphida = (wIbot' * U[1:stc,:,5]) * Wa
    dphids = wDbot' * U[1:stc,:,5]
    dsdx = -dphida ./ dphids
    dsdz = g ./ dphids
    RHS = -rho * g .* NzBot[1,:]' - dPda .* NxBot[1,:]'
    RHS = RHS ./ (NxBot[1,:]' .* dsdx + NzBot[1,:]' .* dsdz)
    U[1,:,6] = (RHS - wDbot[2:stc]' * U[2:stc,:,6]) / wDbot[1]

    #set pressure on top ghost nodes:
    dPda = (wHtop' * U[end-1:-1:end-stc,:,6]) * Wa
    rho = wHtop' * U[end-1:-1:end-stc,:,4]
    dphida = (wItop' * U[end:-1:end-stc+1,:,5]) * Wa
    dphids = wDtop' * U[end:-1:end-stc+1,:,5]
    dsdx = -dphida ./ dphids
    dsdz = g ./ dphids
    RHS = -rho * g .* NzTop[1,:]' - dPda .* NxTop[1,:]'
    RHS = RHS ./ (NxTop[1,:]' .* dsdx + NzTop[1,:]' .* dsdz)
    U[end,:,6] = (RHS - wDtop[2:stc]' * U[end-1:-1:end-stc+1,:,6]) /
        wDtop[1]

    #extrapolate temperature to bottom and top ghost nodes:
    U[1,:,3] = wEbot' * U[2:stc+1,:,3]
    U[end,:,3] = wEtop' * U[end-1:-1:end-stc,:,3]

    #extrapolate density to bottom and top ghost nodes using EOS:
    U[1,:,4] = (Pbar[1,:]+U[1,:,6]) / Rd ./
        (Tbar[1,:]+U[1,:,3]) - rhoBar[1,:]
    U[end,:,4] = (Pbar[end,:]+U[end,:,6]) / Rd ./
        (Tbar[end,:]+U[end,:,3]) - rhoBar[end,:]

    return U, Pbar, rhoBar, Tbar, drhoBarDz, dTbarDz 

end

###########################################################################

#Initial arrays of zeros for storing Runge-Kutta sub-steps.  If we are
#using RK3, then we need only two arrays, but for RK4 we need 4.  Note
#that dUdt and q1 are two different names for the same array.

rks = 3                           #hard-coded: number of Runge-Kutta stages

dUdt = zeros(Float64, nLev+2, nCol, 6)
q1   = dUdt
q2   = zeros(Float64, nLev+2, nCol, 6)

if rks == 4
    q3 = zeros(Float64, nLev+2, nCol, 6)
    q4 = zeros(Float64, nLev+2, nCol, 6)
end

###########################################################################

#This describes the RHS of the system of ODEs in time that will be solved:

function odefun(t, U, dUdt)

    #Preliminaries:

    U, Pbar, rhoBar, Tbar, drhoBarDz, dTbarDz = setGhostNodes(U)

    rhoInv = 1 ./ (rhoBar + U[:,:,4])
    duda   = Da(U[:,:,1])
    duds   = Ds(U[:,:,1])
    dwda   = Da(U[:,:,2])
    dwds   = Ds(U[:,:,2])
    dphids = Ds(U[:,:,5])
    dPds   = Ds(U[:,:,6])

    dsdx = -Da(U[:,:,5]) ./ dphids
    dsdz = g ./ dphids
    divU = (duda + duds .* dsdx) + (dwds .* dsdz)
    uDotGradS = U[:,:,1] .* dsdx + U[:,:,2] .* dsdz

    if verticallyLagrangian
        sDot = zeros(Float64, nLev+2, nCol)
    else
        sDot = uDotGradS
    end

    #Main part:

    dUdt[2:end-1,:,1] = (-U[:,:,1] .* duda - sDot .* duds -
        rhoInv .* (Da(U[:,:,6]) + dPds .* dsdx))[2:end-1,:,1] +
        HV(U[:,:,1])

    dUdt[2:end-1,:,2] = (-U[:,:,1] .* dwda - sDot .* dwds -
        rhoInv .* (dPds .* dsdz) - U[:,:,4] * g .* rhoInv)[2:end-1,:,1] +
        HV(U[:,:,2])

    dUdt[2:end-1,:,3] = (-U[:,:,1] .* Da(U[:,:,3]) - sDot .* Ds(U[:,:,3]) -
        U[:,:,2] .* dTbarDz -
        Rd/Cv * (Tbar+U[:,:,3]) .* divU)[2:end-1,:,1] +
        HV(U[:,:,3])

    dUdt[2:end-1,:,4] = (-U[:,:,1] .* Da(U[:,:,4]) - sDot .* Ds(U[:,:,4]) -
        U[:,:,2] .* drhoBarDz -
        (rhoBar+U[:,:,4]) .* divU)[2:end-1,:,1] +
        HV(U[:,:,4])

    dUdt[2:end-1,:,5] = ((uDotGradS - sDot) .* dphids)[2:end-1,:,1] +
        HV(U[:,:,5] - phiBar)

    return dUdt

end

###########################################################################

#Main time-stepping loop:

# U = setGhostNodes(U)[1]
# figure()
# contourf(xx, zz, U[:,:,6], 20)
# colorbar()
# show()

et = time()

for i in 0:nTimesteps

    global t, U, et, tmp

    if verticallyLagrangian & !contourFromSaved & mod(i,4)==0
        U = setGhostNodes(U)[1]
        U[:,:,1:4] = verticalRemap(U[:,:,1:4], U[:,:,5], phiBar, V)
        U[:,:,5] = phiBar
    end

    if mod(i, Int(round(saveDel/dt))) == 0

        @printf("t = %5d | et=%6.2f | maxAbsRho = %.2e\n",
            Int(round(t)), time()-et, maximum(abs.(U[:,:,4])))
    
        et = time()

        if contourFromSaved
            tmp = readdlm(@sprintf("%s%04d.dat",saveString,Int(round(t))))
            U[:,:,1:4] = reshape(tmp, nLev+2, nCol, 4)
        end 

        if saveArrays | saveContours
            U = setGhostNodes(U)[1]
        end

        if saveArrays
            writedlm(@sprintf("%s%04d.dat",saveString,Int(round(t))),
                U[:,:,1:4][:])
        end

        if saveContours
            contourSomething(U, t)
        end

    end

    if contourFromSaved
        t = t + dt
    else
        if rks == 3
            t, U = rk3!(t, U, odefun, dt, q1, q2)
        elseif rks == 4
            t, U = rk4!(t, U, odefun, dt, q1, q2, q3, q4)
        else
            error("Please use RK3 or RK4 for this problem.")
        end
    end

end

###########################################################################
