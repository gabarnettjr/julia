# Function to compute a Riemann sum.
# Greg Barnett
# Halloween 2020

# f      : The function f(x)
# a      : Left endpoint of the interval
# b      : Right endpoint of the interval
# n      : Number of subintervals to use  ; default = 64
# method : "left", "right", or "midpoint" ; default = "midpoint"

function riemannSum(f, a, b; n=64, method="midpoint")
    dx = (b - a) / n                                # Width of each rectangle
    I = 0                                           # Initial value of the sum
    for j = 0 : n-1                                 # for-loop to build the sum
        if method == "right"
            x = a + (j+1) * dx           # Evaluation point for right-hand rule
        elseif method == "left"
            x = a + j * dx               # Evaluation point for left-hand rule
        elseif method == "midpoint"
            x = a + (j+1/2) * dx         # Evaluation point for midpoint rule
        else
            error("Invalid method flag."
            , "  Please choose \"left\", \"right\", or \"midpoint\".");
        end
        I = I + f(x) * dx            # Add the area of the rectangle to the sum
    end
    return I
end

# HOW TO USE THE FUNCTION IN JULIA:

# From the directory where this file is saved, open the julia prompt,
# and include the file by typing:
# include("riemannSums.jl");

# If you do not know how to use commands to get to the directory where the
# file is saved, then just copy the function from above and then paste it to
# the julia prompt by right-clicking somewhere.

# Once the function is available in julia, you can use it to calculate sums.

# Let's say you want to estimate the area under the function f(x)=x^3 from
# x=1 to x=2, using 16 subintervals and function values at left endpoints.
# You could do it like this:

# f(x) = x^3;  riemannSum(f, 1, 2, n=16, method="left")







