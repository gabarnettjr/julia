# Calculating a few Riemann sums to remember how Julia works.

f(x) = 3*x + 2;
a = 3
b = 5
n = 16
LRM = "M"

function riemannSum(f, a, b, n, LRM)
    dx = (b - a) / n
    I = 0
    for j = 0 : n-1
        if LRM == "R"
            x = a + (j+1) * dx
        elseif LRM == "L"
            x = a + j * dx
        elseif LRM == "M"
            x = a + j * dx + dx/2
        else
            error("Invalid flag LRM.");
        end
        I = I + f(x) * dx
    end
    return I
end

println(riemannSum(f, a, b, n, LRM))

