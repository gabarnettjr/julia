
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os

#####################################################################

# USER INPUT

# Decide whether to plot the eigenvalues
eigenvalues = True

# If the user passes in an argument, assume it's the variable to plot
if len(sys.argv) > 1:
    whatToPlot = sys.argv[1]
else:
    whatToPlot = "rho"

# Directory where the *.txt files are located (in string)
instr = './results/'

# Get the variable string
varstr = instr + whatToPlot

# Directory where the *.png files will end up (out string)
outstr = './figures/'

# Set the contour levels
if whatToPlot == "rho":
    clevels = np.linspace(-.5,.5,41)
elif (whatToPlot == "u") or (whatToPlot == "v"):
    clevels = np.linspace(-1/3, 1/3, 41)

#####################################################################

# Remove old figures if there are any
tmp = os.path.join(os.getcwd(), outstr)
ell = os.listdir(tmp)
for item in ell:
    if item.endswith(".png") :
        os.remove(os.path.join(tmp, item))

# Load the x-coordinates of the points
x = np.array([])
with open(instr + 'x.txt') as f:
    for line in f:
        x = np.hstack((x, np.float(line)))

# Load the y-coordinates of the points
y = np.array([])
with open(instr + 'y.txt') as f:
    for line in f:
        y = np.hstack((y, np.float(line)))

# plot the nodes and save them
fig = plt.figure(figsize = (12, 10))
plt.plot(x, y, '.')
plt.axis('equal')
fig.savefig(outstr + 'nodes.png', bbox_inches = 'tight')


# Load the eigenvalues and plot them, if requested
if eigenvalues:
    e_real = np.array([])
    with open(instr + 'e_real.txt') as f:
        for line in f:
            e_real = np.hstack((e_real, np.float(line)))
    e_imag = np.array([])
    with open(instr + 'e_imag.txt') as f:
        for line in f:
            e_imag = np.hstack((e_imag, np.float(line)))
    fig = plt.figure(figsize = (12, 10))
    plt.plot(e_real, e_imag, '.')
    plt.title('maxReal = {0:.10f}'.format(max(e_real)))
    fig.savefig(outstr + 'eigenvalues.png', bbox_inches = 'tight')

# Get the triangular mesh for plotting the contours
triang = mtri.Triangulation(x, y)

# Initialize the array to store the variable you want to plot
var = np.zeros(np.shape(x))

# Initialize the frame number
frame = 0

# The main loop that plots and saves figures

while True:

    try:

        i = 0
        with open(varstr + "_{0:04d}.txt".format(frame)) as f:
            for line in f:
                var[i] = np.float(line)
                i = i + 1

    except:

        break

    else:
        
        fig = plt.figure(figsize = (12, 10))
        ax = fig.add_subplot(111)
        cs = ax.tricontourf(triang, var, levels = clevels)
        fig.colorbar(cs)
        plt.axis('equal')
        plt.axis([-1,1,-1,1])
        
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_trisurf(tri, rho)
        # ax.set_xlim3d(-1,1)
        # ax.set_ylim3d(-1,1)
        # ax.set_zlim3d(-1,1)
        # plt.show()
        
        fig.savefig(outstr + '{0:04d}'.format(frame) + '.png',
                bbox_inches = 'tight')
        
        plt.close()

        frame = frame + 1







