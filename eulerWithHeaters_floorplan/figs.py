
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os

###############################################################################

# INSTRUCTIONS

# This will run with no command-line inputs.

# HOWEVER, if the user passes in 3 (or 4) arguments, then it assumes they are:
# (1) whether or not to plot vectors (0 or 1)
# (2) whether or not to plot moving points (0 or 1)
# (3) what to plot (vmag, u, v, rho, T, or p)
# (4) (not required) directory to pull files from (default is 'results\')

###############################################################################

# Process the command-line input
if len(sys.argv) >= 4:
    if sys.argv[1] == "1":
        plotVectors = True
    else:
        plotVectors = False
    if sys.argv[2] == "1":
        plotTracers = True
    else:
        plotTracers = False
    whatToPlot = sys.argv[3]
else:
    plotVectors = False
    plotTracers = True
    whatToPlot = "vmag"

# Make it simple by always plotting vmag on same colorbar,
# but using automatic colorbar range for other variables.
if whatToPlot == "vmag":
    autoLevels = False
    clevels = np.arange(0, 22, 2)
else:
    autoLevels = True
    clevels = 9

# Directory where the *.txt files are located (in string)
if len(sys.argv) == 5:
    instr = sys.argv[4]
else:
    instr = 'results/'

# Get the variable string
varstr = instr + whatToPlot

# Directory where the *.png files will end up (out string)
outstr = 'figures/'

###############################################################################

def loadScalar(fileName):
    y = np.array([], float)
    with open(fileName) as f:
        for line in f:
            y = np.hstack((y, float(line)))
    return y[0]

###############################################################################

def loadVector(fileName, isIndex = False, zeroVec = np.zeros((1))):
    if len(zeroVec) > 1:
        y = zeroVec
        i = 0
        with open(fileName) as f:
            for line in f:
                y[i] = float(line)
                i = i + 1
    else:
        if isIndex:
            y = np.array([], int)
        else:
            y = np.array([], float)
        with open(fileName) as f:
            for line in f:
                if isIndex:
                    y = np.hstack((y, int(line) - 1))
                else:
                    y = np.hstack((y, float(line)))
    return y

###############################################################################

# Load the constants and background states

Cp = loadScalar(instr + 'Cp.txt')
Cv = loadScalar(instr + 'Cv.txt')
R  = loadScalar(instr + 'R.txt')

p_0   = loadVector(instr + 'p_0.txt')
rho_0 = loadVector(instr + 'rho_0.txt')
T_0   = loadVector(instr + 'T_0.txt')

###############################################################################

if autoLevels:

    minP   = loadScalar(instr + 'minP.txt')
    maxP   = loadScalar(instr + 'maxP.txt')
    
    minRho = loadScalar(instr + 'minRho.txt')
    maxRho = loadScalar(instr + 'maxRho.txt')
    
    minU   = loadScalar(instr + 'minU.txt')
    maxU   = loadScalar(instr + 'maxU.txt')
    
    minV   = loadScalar(instr + 'minV.txt')
    maxV   = loadScalar(instr + 'maxV.txt')
    
    minT   = loadScalar(instr + 'minT.txt')
    maxT   = loadScalar(instr + 'maxT.txt')

    # Load the maximum velocity magnitude over the whole time interval
    maxVel = loadScalar(instr + 'maxVel.txt')

    # Set the contour levels
    if whatToPlot == "p":
        clevels = np.linspace(minP, maxP, clevels)
    elif whatToPlot == "rho":
        clevels = np.linspace(minRho, maxRho, clevels)
    elif whatToPlot == "u":
        clevels = np.linspace(minU, maxU, clevels)
    elif whatToPlot == "v":
        clevels = np.linspace(minV, maxV, clevels)
    elif whatToPlot == "T":
        clevels = np.linspace(minT, maxT, clevels)
    elif whatToPlot == "vmag":
        clevels = np.linspace(0, maxVel, clevels)
    else:
        sys.exit("Invalid thing to plot.")

###############################################################################

# Remove old figures if there are any
tmp = os.path.join(os.getcwd(), outstr)
ell = os.listdir(tmp)
for item in ell:
    if item.endswith(".png") :
        os.remove(os.path.join(tmp, item))

# Load the x-coordinates of the points
x = loadVector(instr + 'x.txt')

# Load the y-coordinates of the points
y = loadVector(instr + 'y.txt')

# Load the x-coordinates of the NON points
xNon = loadVector(instr + 'xNon.txt')

# Load the y-coordinates of the NON points
yNon = loadVector(instr + 'yNon.txt')

# Load the times that will be plotted:
tSaves = loadVector(instr + 'tSaves.txt')

# Load the approximate node spacing
dr = loadScalar(instr + 'dr.txt')

###############################################################################

# plot the nodes and save the plot

# index of the boundary nodes
bb = loadVector(instr + 'bb.txt', isIndex = True)

# index of the "fan" nodes
ff2 = loadVector(instr + 'ff2.txt', isIndex = True)
try:
    ff3 = loadVector(instr + 'ff3.txt', isIndex = True)
    ff4 = loadVector(instr + 'ff4.txt', isIndex = True)
    ff5 = loadVector(instr + 'ff5.txt', isIndex = True)
except:
    ff3 = []
    ff4 = []
    ff5 = []

fig = plt.figure(figsize = (12, 10))
plt.plot(x, y, 'k.', \
         x[bb], y[bb], 'ks', \
         x[ff2], y[ff2], 'rs', \
         x[ff3], y[ff3], 'rs', \
         x[ff4], y[ff4], 'rs', \
         x[ff5], y[ff5], 'rs', \
         xNon, yNon, 'kx', \
         markersize = 4, fillstyle = 'none')

plt.axis('equal')
# plt.xticks(np.arange(1, 21, 1))
# plt.yticks(np.arange(1, 21, 1))
# plt.grid(True)
fig.savefig(outstr + 'nodes.png', bbox_inches = 'tight')
plt.close()

###############################################################################

# Get the triangular mesh for plotting the contours
triang = mtri.Triangulation(x, y)

# Initialize the array(s) to store the variable(s) you want to plot
var = np.zeros(np.shape(x))
if (plotVectors or (whatToPlot == "vmag")):
    u = np.zeros(np.shape(x))
    v = np.zeros(np.shape(x))
if whatToPlot == "p":
    rho = np.zeros(np.shape(x))
    T = np.zeros(np.shape(x))
if plotTracers:
    xt = np.zeros(np.shape(x))
    yt = np.zeros(np.shape(x))

# Initialize the frame number
frame = 0

# The main loop that plots and saves figures

while True:

    try:

        if plotVectors or (whatToPlot == "vmag"):
            u = loadVector(instr + "u_{0:04d}.txt".format(frame), zeroVec = u)
            v = loadVector(instr + "v_{0:04d}.txt".format(frame), zeroVec = v)
            if (whatToPlot == "vmag"):
                var = np.sqrt(u**2 + v**2)
            if plotVectors:
                u = u / maxVel * dr
                v = v / maxVel * dr
        if whatToPlot == "p":
            rho = loadVector(instr + "rho_{0:04d}.txt".format(frame),
                             zeroVec = rho)
            T = loadVector(instr + "T_{0:04d}.txt".format(frame), zeroVec = T)
            var = rho * R * T
        if (whatToPlot != "vmag") and (whatToPlot != "p"):
            var = loadVector(varstr + "_{0:04d}.txt".format(frame),
                             zeroVec = var)
        if plotTracers:
            xt = loadVector(instr + "xt_{0:04d}.txt".format(frame),
                            zeroVec = xt)
            yt = loadVector(instr + "yt_{0:04d}.txt".format(frame),
                            zeroVec = yt)

    except:

        break

    else:
        
        # fig = plt.figure(figsize = (13, 10))
        fig = plt.figure(figsize = (24, 10))
        # fig = plt.figure(figsize = (18,10))
        # fig = plt.figure()
        ax = fig.add_subplot(111)
        if whatToPlot == "p":
            cs = ax.tricontourf(triang, var-p_0, levels = clevels)
            title = "pressure perturbation"
        elif whatToPlot == "rho":
            cs = ax.tricontourf(triang, var-rho_0, levels = clevels)
            title = "density perturbation"
        elif whatToPlot == "u":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "velocity west-to-east"
        elif whatToPlot == "v":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "velocity south-to-north"
        elif whatToPlot == "T":
            cs = ax.tricontourf(triang, var-T_0, levels = clevels)
            title = "temperature perturbation"
        elif whatToPlot == "vmag":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "air speed (m/s)"
        else:
            sys.exit("Not a valid quantity to  plot.")
        fig.colorbar(cs)

        plt.axis('equal')
        plt.axis(np.array([np.min(x), np.max(x), np.min(y), np.max(y)]) +
                 np.array([-1/5,      1/5,       -1/5,      1/5]))

        # plt.plot(x[bb], y[bb], 'k.', markersize = 5)

        plt.plot(x[ff2], y[ff2], 'r.', markersize = 5)
        plt.plot(x[ff3], y[ff3], 'r.', markersize = 5)
        plt.plot(x[ff4], y[ff4], 'r.', markersize = 5)
        plt.plot(x[ff5], y[ff5], 'r.', markersize = 5)

        # k  =  5,  7,  9, 11, 13, 15
        # ms = 20, 14, 12,  9,  8,  7
        plt.plot(xNon, yNon, 'ks', markersize = 6)
        
        if plotVectors:
            plt.plot(np.vstack((x-u/2, x+u/2)),
                     np.vstack((y-v/2, y+v/2)), 'k')
        if plotTracers:
            plt.plot(xt, yt, 'k.', markersize = 1)

        plt.title(title + ',    t = {0:.5f} seconds'.format(tSaves[frame]))
        plt.xlabel('meters')
        
        fig.savefig(outstr + '{0:04d}'.format(frame) + '.png',
                    bbox_inches = 'tight')
        
        plt.close()

        frame = frame + 1

