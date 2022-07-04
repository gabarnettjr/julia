
function heatPrime(t::Float64; td::Float64 = 20., a::Float64 = 10000.)

    return 3 * td * a * t ^ 2 / (1 + a * t ^ 3) ^ 2

    # return -16 * td * (4*(t-1/2))^3 * exp(-(4*(t-1/2))^4)

end

