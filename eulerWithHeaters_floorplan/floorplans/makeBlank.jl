
# This makes a blank slate of mostly ones (where air is allowed to flow), but with an outer border of zeros (solid walls).  It creates a good blank slate for creating your own floorplan.  To modify it, just find where you want to put solid walls and change the ones to zeros there.
# Fans that blow west-to-east are indicated with the number 2.
# Fans that blow south-to-north are indicated with the number 3.
# Fans that blow east-to-west are indicated with the number 4.
# Fans that blow north-to-south are indicated with the number 5.

# Call this from the command prompt using three inputs.
# The first input is the desired name of the file (end it with .txt).
# The second input is the desired number of rows.
# The third input is the desired number of columns.

using DelimitedFiles

fileName = ARGS[1]
len = parse(Int, ARGS[2])
wid = parse(Int, ARGS[3])

M = ones(Int, len, wid)

M[[1,end], :] .= 0
M[:, [1,end]] .= 0

# M[[4,end-3], 4:end-3] .= 0
# M[4:end-3, [4,end-3]] .= 0

io = open(fileName, "w")
writedlm(io, M, ' ')
close(io)
