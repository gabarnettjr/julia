using Printf
using PyPlot

include("rk.jl")
include("phs1.jl")

###########################################################################

testCase = "risingBubble"

verticalCoordinate = "pressure"
verticallyLagrangian = true

#Choose 0, 1, 2, 3, or 4:
refinementLevel = 1

saveArrays       = true
saveContours     = true
contourFromSaved = false
plotNodesAndExit = false

#Choose "u", "w", "T", "rho", "phi", "P", "theta", "pi", or "phi":
whatToPlot = "phi"

#Choose either a number of contours, or a range of contours:
contours = 20

###########################################################################

#Get string for saving results:

saveString = string("./nonhydroXZresults/", testCase, "_",
    verticalCoordinate, "_")

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
    tmp = walkdir("nonhydroXZresults")
    for (roots, dirs, files) in tmp
        for file in files
            if endswith(file, ".png") | endswith(file, ".PNG")
                if Sys.islinux()
                    # print(string("nonhydroXZresults","/",file,"\n"))
                    rm(string("nonhydroXZresults","/",file))
                elseif Sys.iswindows()
                    # print(string("nonhydroXZresults","\\",file,"\n"))
                    rm(string("nonhydroXZresults","\\",file))
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

    function inverseExnerPressure(pi)
        return (1 - pi) * Cp*th0/g
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
        return 2 * exp.(-1e-6*(x.^2 + (z-4000).^2))
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

    if verticalCoordinate == "height"

        ds = (top - bottom) / nLev
        s = range(bottom-ds/2, stop=top+ds/2, length=nLev+2)

        pSurf = 0
        sTop = 0
        pTop = 0

    else

        piSurf = exnerPressure(zSurf)
        piTop  = exnerPressure(top)
        pSurf = Po * piSurf .^ (Cp/Rd)     #hydrostatic pressure at surface
        pTop  = Po * piTop  ^ (Cp/Rd)          #hydrostatic pressure at top
        sTop = pTop / Po                      #value of s on upper boundary
        ds = (1 - sTop) / nLev
        s = range(sTop-ds/2, stop=1+ds/2, length=nLev+2)

    end

    return s, ds, sTop, pTop, pSurf

end

###########################################################################

s, ds, sTop, pTop, pSurf = getSvalues()

xx = repeat(transpose(x), nLev+2, 1)
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

Wa = transpose(Wa)
Whva = transpose(Whva)

phs = 5                                              #vertical PHS exponent
pol = 3                        #highest degree polynomial in vertical basis
stc = 7                                              #vertical stencil size
if verticalCoordinate == "pressure"
    alp = -2^-20 * 300                    #vertical dissipation coefficient
else
    alp = -2^-10 * 300               #much larger in height coordinate case
end
Ws   = getDM(z=s, x=s,          m=1,     phs=phs, pol=pol, stc=stc)
Whvs = getDM(z=s, x=s[2:end-1], m=pol+1, phs=phs, pol=pol, stc=stc)
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
    return U[2:end-1,] * Whva + Whvs * U
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
