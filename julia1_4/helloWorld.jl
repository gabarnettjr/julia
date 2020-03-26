println("Hello, World!")

println("That was not so easy.")

println("What else can you do with this?")

function addThreeNumbers( x, y, z )
    w = x + y + z
    return w
end

function subtractThreeNumbers( x, y, z )
    w = x - y - z
    s = x - y + z
    return ( w, s )
end

println(addThreeNumbers(1,2,3))

println(subtractThreeNumbers(6,5,4))

