
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os

#####################################################################

# USER INPUT

# If the user passes in an argument, assume it's the variable to plot
if len(sys.argv) > 1:
    whatToPlot = sys.argv[1]
else:
    whatToPlot = "rho"

# Directory where the *.txt files are located (in string)
instr = 'results/'

# Get the variable string
varstr = instr + whatToPlot

# Directory where the *.png files will end up (out string)
outstr = 'figures/'

# Set the contour levels
clevels = 20
# if whatToPlot == "rho":
#     clevels = np.arange(-105, 65, 10) * 1e-5
# elif whatToPlot == "u":
#     clevels = np.arange(-12.5, 13.5, 1)
# elif whatToPlot == "v":
#     clevels = np.arange(-10.5, 10.5, 1)
# elif whatToPlot == "e":
#     clevels = np.arange(-75, 45, 10)
# else:
#     sys.exit("Invalid thing to plot.")

#####################################################################

# Remove old figures if there are any
tmp = os.path.join(os.getcwd(), outstr)
ell = os.listdir(tmp)
for item in ell:
    if item.endswith(".png") :
        os.remove(os.path.join(tmp, item))

# Load the x-coordinates of the points
x = np.array([], float)
with open(instr + 'x.txt') as f:
    for line in f:
        x = np.hstack((x, np.float(line)))

# Load the y-coordinates of the points
y = np.array([], float)
with open(instr + 'y.txt') as f:
    for line in f:
        y = np.hstack((y, np.float(line)))

# Load the radius of the boundary of the disk
radius = np.array([], float)
with open(instr + 'radius.txt') as f:
    for line in f:
        radius = np.hstack((radius, np.float(line)))
radius = radius[0]

###############################################################################

# plot the nodes and save the plot

ind_ghost = np.array([], int)
with open(instr + 'ind_ghost.txt') as f:
    for line in f:
        ind_ghost = np.hstack((ind_ghost, np.int(line) - 1))

ind_noGhost = np.array([], int)
with open(instr + 'ind_noGhost.txt') as f:
    for line in f:
        ind_noGhost = np.hstack((ind_noGhost, np.int(line) - 1))

ind_outer = np.array([], int)
with open(instr + 'ind_outer.txt') as f:
    for line in f:
        ind_outer = np.hstack((ind_outer, np.int(line) - 1))

ind_fan = np.array([], int)
with open(instr + 'ind_fan.txt') as f:
    for line in f:
        ind_fan = np.hstack((ind_fan, np.int(line) - 1))

th = np.linspace(0, 2*np.pi, 250)

fig = plt.figure(figsize = (12, 10))
plt.plot(x, y, 'k.', \
         radius*np.cos(th), radius*np.sin(th), 'k', \
         # x[indA_inflow], y[indA_inflow], 'ro', \
         # x[indB_inflow], y[indB_inflow], 'go', \
         # x[indC_inflow], y[indC_inflow], 'yo', \
         # x[indA_outflow], y[indA_outflow], 'rd', \
         # x[indB_outflow], y[indB_outflow], 'gd', \
         # x[indC_outflow], y[indC_outflow], 'yd', \
         x[ind_outer], y[ind_outer], 'gs', \
         x[ind_ghost], y[ind_ghost], 'ys', \
         x[ind_fan], y[ind_fan], 'bs', \
         markersize = 5, fillstyle = 'none')
plt.axis('equal')
fig.savefig(outstr + 'nodes.png', bbox_inches = 'tight')
plt.close()

###############################################################################

# Load the background states

rho_0 = np.array([], float)
with open(instr + 'rho_0.txt') as f:
    for line in f:
        rho_0 = np.hstack((rho_0, np.float(line)))

u_0 = np.array([], float)
with open(instr + 'u_0.txt') as f:
    for line in f:
        u_0 = np.hstack((u_0, np.float(line)))

v_0 = np.array([], float)
with open(instr + 'v_0.txt') as f:
    for line in f:
        v_0 = np.hstack((v_0, np.float(line)))

e_0 = np.array([], float)
with open(instr + 'e_0.txt') as f:
    for line in f:
        e_0 = np.hstack((e_0, np.float(line)))

###############################################################################

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
        if whatToPlot == "rho":
            cs = ax.tricontourf(triang, var-rho_0, levels = clevels)
        elif whatToPlot == "u":
            cs = ax.tricontourf(triang, var-u_0, levels = clevels)
        elif whatToPlot == "v":
            cs = ax.tricontourf(triang, var-v_0, levels = clevels)
        elif whatToPlot == "e":
            cs = ax.tricontourf(triang, var-e_0, levels = clevels)
        else:
            sys.exit("Not a valid thing to  plot.")
        fig.colorbar(cs)
        plt.axis('equal')
        plt.axis(12 * np.array([-1,1,-1,1]))

        plt.plot(x[ind_fan], y[ind_fan], 'k--')
        plt.plot(radius*np.cos(th), radius*np.sin(th), 'k')
        
        fig.savefig(outstr + '{0:04d}'.format(frame) + '.png',
                bbox_inches = 'tight')
        
        plt.close()

        frame = frame + 1







