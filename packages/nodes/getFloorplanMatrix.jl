
# Define the layout of the floor using zeros to mark the boundaries, and using
# ones to mark the space where air is allowed to flow.  Also, use the number 2
# to indicate where a fan will be forcing the air left-to-right.

using DelimitedFiles

function getFloorplanMatrix(; plan = "square", useTxt = true)

    if useTxt
        M = readdlm("floorplans/" * plan * ".txt", Int);
        M = reverse(M, dims=1)
        (len, wid) = size(M);
        return M, len, wid
    end

    if plan == "square" || plan == "midSquare" || plan == "midSquare13" ||
    plan == "isoSquare"
        if plan == "midSquare13"
            len = 13
        else
            len = 9
        end
        wid = 18
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        if plan == "isoSquare"
            M[2:end-1,3] .= 2
        else
            M[4:end-3, 3] .= 2
            M[end-2, 3:end-2] .= 0
            M[3, 3:end-2] .= 0
        end
        if plan == "square" || plan == "isoSquare"
            M[Int((len+1)/2), 5] = 0
        else
            M[Int((len+1)/2), Int(wid/2)] = 0
        end
        M = reverse(M, dims=1)
    elseif plan == "default" || plan == "linz"
        len = 12
        wid = 16
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        M[4:9, 4] = [0; 2; 0; 0; 0; 0]
        # M[5, 4] = 1 # replace fan with empty space
        M[[4,9], [5,6,11,12]] .= 0
        M[[4,5,6,9], [7,10]] .= 0
        M[[6,9], [8,9]] .= 0
        M[4:9, 13] = [0; 0; 0; 1; 0; 0]
        if plan == "linz"
            M[7, 13] = 0
        end
        M[5,3] = 3
        M = reverse(M, dims=1)
    elseif plan == "noWalls" || plan == "noWallsSymmetric" || plan == "tunnel"
        len = 8
        wid = 8
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        if plan == "noWalls"
            M[5:6, 4] .= 2
            M[5:6, 2:3] .= 3
        elseif plan == "noWallsSymmetric"
            M[4:5, 4] .= 2
            M[4:5, 2:3] .= 3
        end
        if plan == "tunnel"
            M[[3,6], 3:6] .= 0
            M[4:5,3] .= 2
        end
    elseif plan == "I" || plan == "I2"
        len = 20
        wid = 20
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        M[6:7, 3] .= 3
        M[6:7, 4] .= 2
        M[5, 4] = 0
        M[4, 4:17] .= 0
        M[6:9, 17] .= 0
        M[9, 12:16] .= 0
        M[10:12, 12] .= 0
        M[12, 13:17] .= 0
        M[13, 17] = 0
        M[16, 17] = 0
        M[17, 4:17] .= 0
        M[12:15, 4] .= 0
        M[12, 5:9] .= 0
        M[9:11, 9] .= 0
        M[9, 4:8] .= 0
        M[8, 4] = 0
        if plan == "I2"
            M[5,17] = 1
        end
        M = reverse(M, dims=1)
    elseif plan == "vinny"
        len = 28
        wid = 33
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        M[4:end-3, 30] .= 0
        M[4, 4:end-3] .= 0
        M[4:end-3, 4] .= 0
        M[end-3, 4:end-3] .= 0
        M[end-11, end-8:end-3] .= 0
        M[end-11, end-13:end-11] .= 0
        M[end-10:end-9, end-13] .= 0
        M[end-5:end-3, end-13] .= 0
        M[10, end-7:end-3] .= 0
        M[11:12, end-7] .= 0
        M[15:16, end-7] .= 0
        M[5:10, 19] .= 0
        M[10, 20:23] .= 0
        M[4:10, 14] .= 0
        M[10, 11:14] .= 0
        M[10, 5:8] .= 0
        M[11, 8] = 0
        M[14, 8] = 0
        M[15, 4:8] .= 0
        M[18, 4:9] .= 0
        M[18, 13:14] .= 0
        M[19:24, 14] .= 0
        M[12, 4] = 2
        M[12, 3] = 3
        M[end-3, 8:10] .= 1
        M[end-3, end-8:end-7] .= 1
        M[4, end-8:end-7] .= 1
        M[4, 9:10] .= 1
        M = reverse(M, dims=1)
    elseif plan == "vinny2"
        len = 25
        wid = 25
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        # walls
        M[4:22, 22] .= 0
        M[4, 4:21] .= 0
        M[4:22, 4] .= 0
        M[22, 4:21] .= 0
        M[15, 16:21] .= 0
        M[16:21, 16] .= 0
        M[9:14, 18] .= 0
        M[5:9, 15] .= 0
        M[9, 16:21] .= 0
        M[9, 5:11] .= 0
        M[5:8, 11] .= 0
        M[10:16, 8] .= 0
        M[17, 5:11] .= 0
        M[18:21, 11] .= 0
        # windows
        M[4, 18:19] .= 1
        M[4, 7:8] .= 1
        M[22, 7:8] .= 1
        M[22, 19:20] .= 1
        # doorways
        M[18:19, 16] .= 1
        M[12, 18] = 1
        M[9, 16] = 1
        M[4, 13] = 1
        M[9, 10] = 1
        M[12:13, 8] .= 1
        M[17, 10] = 1
        # fans
        M[12:13, 4] .= 2
        # tracer
        M[12:13, 2:3] .= 3
        M = reverse(M, dims=1)
    elseif plan == "loop" || plan == "loop2"
        len = 24
        wid = 22
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        M[17, 4:6] .= 0
        M[18:21, 6] .= 0
        M[21, 7:19] .= 0
        M[4:20, 19] .= 0
        M[4, 6:18] .= 0
        M[5:8, 6] .= 0
        M[8, 4:5] .= 0
        M[14, 4:9] .= 0
        M[15:18, 9] .= 0
        M[18, 10:16] .= 0
        M[7:17, 16] .= 0
        M[7, 9:15] .= 0
        M[8:11, 9] .= 0
        M[11, 4:8] .= 0
        if plan == "loop"
            M[12:13, 4] .= 0
        end
        M[15:16, 4] .= 2
        M[15:16, 2:3] .= 3
        M = reverse(M, dims=1)
    elseif plan == "vinny3" || plan == "vinny4" || plan == "vinny5"
        len = 25
        wid = 25
        M = zeros(Int, len, wid)
        M[2:end-1, 2:end-1] .= 1
        # outside walls
        M[22, 20:22] .= 0
        M[19:21, 22] .= 0
        M[13:17, 22] .= 0
        M[7:11, 22] .= 0
        M[4:5, 22] .= 0
        M[4, 20:21] .= 0
        M[4, 14:18] .= 0
        M[4, 9:12] .= 0
        M[4, 4:7] .= 0
        M[5, 4] = 0
        M[6, 4] = 2
        M[5:7, 2:3] .= 3
        M[7:10, 4] .= 0
        M[11, 4] = 2
        M[10:12, 2:3] .= 3
        M[12:14, 4] .= 0
        M[15,4] = 2
        M[14:16, 2:3] .= 3
        M[16:19, 4] .= 0
        M[20, 4] = 2
        M[19:21, 2:3] .= 3
        M[21:22, 4] .= 0
        M[22, 4:7] .= 0
        M[22, 9:18] .= 0
        # inside walls
        M[15, 16:21] .= 0
        M[16:18, 16] .= 0
        M[21, 16] = 0
        M[9, 17:21] .= 0
        M[10:11, 18] .= 0
        M[13:14, 18] .= 0
        M[5:9, 15] .= 0
        M[5:9, 11] .= 0
        M[9, 5:9] .= 0
        M[10:11, 9] .= 0
        M[13, 9] = 0
        M[14, 5:9] .= 0
        M[15, 8] = 0
        M[17, 5:11] .= 0
        M[18, 11] = 0
        M[21, 11] = 0
        if plan == "vinny4" || plan == "vinny5"
            # close north and south windows
            M[22,8] = 0
            M[16, 8] = 0
            M[15, 8] = 1
            M[22, 19] = 0
            M[4, 8] = 0
            M[4, 19] = 0
            if plan == "vinny5"
                # close the back door
                M[4, 13] = 0
            end
        end
        M = reverse(M, dims=1)
    else
        error("That is not a valid floorplan.")
    end

    return M, len, wid

end

