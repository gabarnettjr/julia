
function fanPrime(t::Float64; fs::Float64 = 10., a::Float64 = 10000.)

    return 3 * fs * a * t ^ 2 / (1 + a * t ^ 3) ^ 2

end

