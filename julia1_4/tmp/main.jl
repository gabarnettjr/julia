# x = 0;
# 
# function myForLoop(x)
#     for i in 1 : 100
#         x = x + 1;
#     end
#     return x
# end
# 
# x = myForLoop(x)

####################################################

using Primes

a = 9

b = 18

c = 27

function doStuff( a, b, c; N=101, k=0, x=[] )

    for i in range( -(N-1)/2, stop=(N-1)/2, step=1 )
        tmp = c * i + b
        if mod( tmp, a ) == 0
            push!( x, tmp/a )
            k = k + 1
        end
    end

    return x, k

end

x, k = doStuff( a, b, c )

println()
println( k )
println()
println( x )
println()
println( x[2:k] - x[1:k-1] )
println()

