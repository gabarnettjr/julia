using Plots

# The t values where data is known:
t = range(0, stop=40, step=5)

# The known data (function) values at each t value:
b = [ 1; 5; 8; 12; 13; 15; 35; 59; 101 ]

# The matrix A for the over-determined linear system Ax=b:
A = zeros(length(t), 4)
A[:,1] = t .^ 3
A[:,2] = t .^ 2
A[:,3] = t
A[:,4] = ones( size(t) )

# Solve the normal equations using the backslash operator:
x = A \ b

# The least-squares cubic polynomial function:
function P(t; c=x)
    c[1] * t.^3 .+ c[2] * t.^2 .+ c[3] * t .+ c[4] * ones(size(t))
end

# A dense set of t-values between 0 and 40:
T = range(0, stop=40, length=80)

# Plot showing the polynomial and the given data values:
p1 = plot(T, P(T), lw=2, label="cubic")
scatter!(t, b, label="given data", title="Least-Squares Cubic")

# Plot showing the 9 values of the residual vector:
p2 = plot(t, P(t)-b, legend=false)
scatter!(t, P(t)-b, title="Residual")

# Put both plots together in one figure:
plot(p1, p2, layout=(1,2))

