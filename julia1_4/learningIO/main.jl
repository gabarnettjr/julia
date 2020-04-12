#####################################################################

# Old school way without using DelimitedFiles

# io = open("myfile.txt", "w")
# 
# writedlm(io, "Hello, world!");
# 
# close(io);
# 
# io = open("myfile.txt", "r");
# 
# s = read(io, String);
# print(s)
# 
# close(io);
# 
# rm("myfile.txt");

#####################################################################

using DelimitedFiles

#####################################################################

# Remove any old versions that might already be saved

if isfile("x.txt")
    rm("x.txt");
end

if isfile("A.txt")
    rm("A.txt");
end

if isfile("B.txt")
    rm("B.txt");
end

#####################################################################

# Create Chebychev nodes and write them to a text file

N = 16;

theta = range(0, stop = pi, length = N);
x = cos.(theta);
x = reverse(x);

io = open("x.txt", "w");
writedlm(io, x, ' ');
close(io);

#####################################################################

# Create a random matrix and save it to a text file

A = rand(5,3);

io = open("A.txt", "w");
writedlm(io, A, ' ');
close(io);

#####################################################################

# Create a random array and save it to a text file

B = rand(3,2,4);

io = open("B.txt", "w");
writedlm(io, B, ' ');
close(io);

#####################################################################

# Read in the contents of the files and make them the right size

io = open("x.txt", "r");
new_x = readdlm(io, ' ', Float64);
close(io);
new_x = reshape(new_x, size(x))

io = open("A.txt", "r");
new_A = readdlm(io, ' ', Float64);
close(io);
new_A = reshape(new_A, size(A));

io = open("B.txt", "r");
new_B = readdlm(io, ' ', Float64);
close(io);
new_B = reshape(new_B, size(B));

println(maximum(abs.(x - new_x)))
println(maximum(abs.(A - new_A)))
print(maximum(abs.(B - new_B)))

#####################################################################

