using PyPlot

p = -1 .+ 2*rand(Float64,58,2)
x = p[:,1];  y = p[:,2]
figure()
scatter( x, y )
axis(:image)
axis([-1,1,-1,1])
show()

x = range( -1, stop=1, step=.05 )
y = range( -1, stop=1, length=50 )
z = zeros( length(y), length(x) )
for i in 1 : length(y)
    for j in 1 : length(x)
        z[i,j] = cos( pi*x[j] ) .* sin( pi*y[i] )
    end
end
figure()
contourf( x, y, z, range(-1.05,stop=1.05,step=.1) )
colorbar()
axis(:image)
show()
