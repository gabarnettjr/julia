###########################################################################

function constants(testCase)

    Cp = 1004                #specific heat of dry air at constant pressure
    Cv = 717                   #specific heat dry of air at constant volume
    Rd = Cp - Cv                                  #gas constant for dry air
    g  = 9.81                                       #gravitational constant
    Po = 10^5                                #reference pressure in Pascals

    if testCase == "scharMountainWaves"
        th0 = 288
    else
        th0 = 300                #reference potential temperature in kelvin
    end

    N   = .01        #Brunt-Vaisala frequency for inertia gravity wave case

    return Cp, Cv, Rd, g, Po, th0, N

end

###########################################################################

function domainParameters(testCase, refinementLevel, g, Cp, th0)
    """
    Get some test-specific parameters, such as the size of the domain, the
    horizontal node spacing dx, the number of vertical levels nLev, the
    time between saves saveDel, the topogrophy function zSurf(x), and the
    initial perturbation in potential temperature thetaPtb(x,z):
    """
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

        function zSurf(x)
            return 1000 * sin.(2*pi / 10000 * x)
        end

        function thetaPtb(x, z)
            R = 1500
            xc = 0
            zc = 3000
            r = sqrt.((x-xc).^2 + (z-zc).^2)
            ind = r < R
            C0bubble = zeros(size(x))
            C0bubble[ind] = 2 * ( 1 .- r[ind]./R )
            return C0bubble
            # return 2 * exp(-1e-6*(x.^2 + (z-4000).^2))
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
        saveDel = 75

        function zSurf(x)
            return zeros(size(x))
        end

        function thetaPtb(x, z)
            xc = 0
            zc = 3000
            xr = 4000
            zr = 2000
            rTilde = sqrt.(((x - xc) ./ xr).^2 + ((z - zc) ./ zr).^2)
            Tptb = zeros(size(x))
            ind = rTilde <= 1
            Tptb[ind] = -15/2 * (1 .+ cos.(pi * rTilde[ind]))
            piBar = 1 .- g / Cp / th0 * z
            return Tptb / piBar 
        end

    elseif testCase == "inertiaGravityWaves"

        left    = -150000
        right   = 150000
        bottom  = 0
        top     = 10000
        dx      = 1000 / 2^refinementLevel
        nLev    = 10 * 2^refinementLevel
        dt      = 1 / 2^refinementLevel
        tf      = 3000
        saveDel = 250

        function zSurf(x)
            return zeros(size(x))
        end

        function thetaPtb(x, z)
            thetaC = .01
            hC = 10000
            aC = 5000
            xC = -50000
            return thetaC * sin.(pi*z./hC) ./ (1 .+ ((x - xC) ./ aC).^2)
        end

    elseif testCase == "steadyState"

        left    = -5000
        right   = 5000
        bottom  = 0
        top     = 10000
        dx      = 1000 / 2^refinementLevel
        nLev    = 10 * 2^refinementLevel
        dt      = 1 / 2^refinementLevel
        tf      = 5000
        saveDel = 500

        function zSurf(x)
            return 1000 * exp(-1e-6 * x.^2)
            # return 1000 * sin.(2*pi ./ 10000 .* x)
        end

        function thetaPtb(x, z)
            return zeros(size(x))
        end

    elseif testCase == "scharMountainWaves"

        left = -100000
        right = 100000
        bottom = 0
        top = 19500
        dx = 500
        nLev = 65
        dt = 1
        tf = 20000
        saveDel = 2000

        function zSurf(x)
            h0 = 250
            a = 5000
            lam = 4000
            return h0 * exp(-(x./a).^2) .* (cos.(pi*x./lam)).^2
        end

        function thetaPtb(x, z)
            return zeros(size(x))
        end

    else

        error("Invalid testCase string.  Please choose 'risingBubble', ",
           "'densityCurrent', 'inertiaGravityWaves', '', or 'steadyState'.")

    end

    return left, right, bottom, top, dx, nLev, dt, tf, saveDel,
        zSurf, thetaPtb

end

###########################################################################
