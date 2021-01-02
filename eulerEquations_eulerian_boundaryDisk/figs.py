
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os

#####################################################################

# USER INPUT

autoLevels = True
clevels = np.arange(-.025, 1.075, .05)

# It can automatically choose the contour levels if this is True
if autoLevels:
    clevels = 20

# If the user passes in 3 arguments, assume the first argument is whether
# or not to plot vectors (0 or 1), the second argument is whether or not to
# plot moving points (0 or 1), and the third argument is what to plot, which
# can be u, v, w, rho, q, or vmag.
if len(sys.argv) == 4:
    plotVectors = sys.argv[1]
    if plotVectors == "1":
        plotVectors = True
    else:
        plotVectors = False
    plotTracers = sys.argv[2]
    if plotTracers == "1":
        plotTracers = True
    else:
        plotTracers = False
    whatToPlot = sys.argv[3]
else:
    plotVectors = False
    plotTracers = False
    whatToPlot = "u"

# Directory where the *.txt files are located (in string)
instr = 'results/'

# Get the variable string
varstr = instr + whatToPlot

# Directory where the *.png files will end up (out string)
outstr = 'figures/'

#####################################################################

# Load the constants and background states

Cv = np.array([], float)
with open(instr + 'Cv.txt') as f:
    for line in f:
        Cv = np.hstack((Cv, np.float(line)))
Cv = Cv[0]

R = np.array([], float)
with open(instr + 'R.txt') as f:
    for line in f:
        R = np.hstack((R, np.float(line)))
R = R[0]

rho_0 = np.array([], float)
with open(instr + 'rho_0.txt') as f:
    for line in f:
        rho_0 = np.hstack((rho_0, np.float(line)))

e_0 = np.array([], float)
with open(instr + 'e_0.txt') as f:
    for line in f:
        e_0 = np.hstack((e_0, np.float(line)))

p_0 = rho_0 * R * (e_0 / Cv)

###############################################################################

if autoLevels:

    minRho = np.array([], float)
    with open(instr + 'minRho.txt') as f:
        for line in f:
            minRho = np.hstack((minRho, np.float(line)))
    minRho = minRho[0]
    
    maxRho = np.array([], float)
    with open(instr + 'maxRho.txt') as f:
        for line in f:
            maxRho = np.hstack((maxRho, np.float(line)))
    maxRho = maxRho[0]
    
    minU = np.array([], float)
    with open(instr + 'minU.txt') as f:
        for line in f:
            minU = np.hstack((minU, np.float(line)))
    minU = minU[0]
    
    maxU = np.array([], float)
    with open(instr + 'maxU.txt') as f:
        for line in f:
            maxU = np.hstack((maxU, np.float(line)))
    maxU = maxU[0]
    
    minV = np.array([], float)
    with open(instr + 'minV.txt') as f:
        for line in f:
            minV = np.hstack((minV, np.float(line)))
    minV = minV[0]
    
    maxV = np.array([], float)
    with open(instr + 'maxV.txt') as f:
        for line in f:
            maxV = np.hstack((maxV, np.float(line)))
    maxV = maxV[0]
    
    minE = np.array([], float)
    with open(instr + 'minE.txt') as f:
        for line in f:
            minE = np.hstack((minE, np.float(line)))
    minE = minE[0]
    
    maxE = np.array([], float)
    with open(instr + 'maxE.txt') as f:
        for line in f:
            maxE = np.hstack((maxE, np.float(line)))
    maxE = maxE[0]

    minQ = np.array([], float)
    with open(instr + 'minQ.txt') as f:
        for line in f:
            minQ = np.hstack((minQ, np.float(line)))
    minQ = minQ[0]
    
    maxQ = np.array([], float)
    with open(instr + 'maxQ.txt') as f:
        for line in f:
            maxQ = np.hstack((maxQ, np.float(line)))
    maxQ = maxQ[0]

    minP = np.array([], float)
    with open(instr + 'minP.txt') as f:
        for line in f:
            minP = np.hstack((minP, np.float(line)))
    minP = minP[0]
    
    maxP = np.array([], float)
    with open(instr + 'maxP.txt') as f:
        for line in f:
            maxP = np.hstack((maxP, np.float(line)))
    maxP = maxP[0]

    # Load the maximum velocity magnitude over the whole time interval
    maxVel = np.array([], float)
    with open(instr + 'maxVel.txt') as f:
        for line in f:
            maxVel = np.hstack((maxVel, np.float(line)))
    maxVel = maxVel[0]

    # Set the contour levels
    if whatToPlot == "rho":
        clevels = np.linspace(minRho, maxRho, clevels)
    elif whatToPlot == "u":
        clevels = np.linspace(minU, maxU, clevels)
    elif whatToPlot == "v":
        clevels = np.linspace(minV, maxV, clevels)
    elif whatToPlot == "e":
        clevels = np.linspace(minE, maxE, clevels)
    elif whatToPlot == "q":
        clevels = np.linspace(minQ, maxQ, clevels)
    elif whatToPlot == "vmag":
        clevels = np.linspace(0, maxVel, clevels)
    elif whatToPlot == "p":
        clevels = np.linspace(minP, maxP, clevels)
    else:
        sys.exit("Invalid thing to plot.")

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

# Load the x-coordinates of the NON points
xNon = np.array([], float)
with open(instr + 'xNon.txt') as f:
    for line in f:
        xNon = np.hstack((xNon, np.float(line)))

# Load the y-coordinates of the NON points
yNon = np.array([], float)
with open(instr + 'yNon.txt') as f:
    for line in f:
        yNon = np.hstack((yNon, np.float(line)))

# Load the times that will be plotted:
tSaves = np.array([], float)
with open(instr + 'tSaves.txt') as f:
    for line in f:
        tSaves = np.hstack((tSaves, np.float(line)))

# Load the radius of the disk domain
radius = np.array([], float)
with open(instr + 'radius.txt') as f:
    for line in f:
        radius = np.hstack((radius, np.float(line)))
radius = radius[0]

# Load the approximate node spacing
dr = np.array([], float)
with open(instr + 'dr.txt') as f:
    for line in f:
        dr = np.hstack((dr, np.float(line)))
dr = dr[0]

###############################################################################

# plot the nodes and save the plot

# index of the boundary nodes
bb = np.array([], int)
with open(instr + 'bb.txt') as f:
    for line in f:
        bb = np.hstack((bb, np.int(line) - 1))

# index of the "fan" nodes
ff = np.array([], int)
with open(instr + 'ff.txt') as f:
    for line in f:
        ff = np.hstack((ff, np.int(line) - 1))

th = np.linspace(0, 2*np.pi, 250)

fig = plt.figure(figsize = (12, 10))
plt.plot(x, y, 'k.', \
         radius*np.cos(th), radius*np.sin(th), 'k', \
         x[bb], y[bb], 'ks', \
         x[ff], y[ff], 'rs', \
         xNon, yNon, 'kx', \
         markersize = 4, fillstyle = 'none')
plt.axis('equal')
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
if plotTracers:
    xt = np.zeros(np.shape(x))
    yt = np.zeros(np.shape(x))
if (whatToPlot == "p"):
    rho = np.zeros(np.shape(x))
    e = np.zeros(np.shape(x))

# Initialize the frame number
frame = 0

# The main loop that plots and saves figures

while True:

    try:

        if (plotVectors or (whatToPlot == "vmag")):
            i = 0
            with open("results/u" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    u[i] = np.float(line)
                    i = i + 1
            i = 0
            with open("results/v" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    v[i] = np.float(line)
                    i = i + 1
            if (whatToPlot == "vmag"):
                var = np.sqrt(u**2 + v**2)
            if plotVectors:
                u = u / maxVel * dr
                v = v / maxVel * dr
        elif (whatToPlot == "p"):
            i = 0
            with open("results/rho" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    rho[i] = np.float(line)
                    i = i + 1
            i = 0
            with open("results/e" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    e[i] = np.float(line)
                    i = i + 1
            var = rho * R * (e / Cv)
        else:
            i = 0
            with open(varstr + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    var[i] = np.float(line)
                    i = i + 1
        if plotTracers:
            i = 0
            with open("results/xt" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    xt[i] = np.float(line)
                    i = i + 1
            i = 0
            with open("results/yt" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    yt[i] = np.float(line)
                    i = i + 1

    except:

        break

    else:
        
        fig = plt.figure(figsize = (12, 10))
        ax = fig.add_subplot(111)
        if whatToPlot == "rho":
            cs = ax.tricontourf(triang, var-rho_0, levels = clevels)
            title = "density perturbation"
        elif whatToPlot == "u":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "horizontal velocity"
        elif whatToPlot == "v":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "vertical velocity"
        elif whatToPlot == "e":
            cs = ax.tricontourf(triang, var-e_0, levels = clevels)
            title = "energy perturbation"
        elif whatToPlot == "q":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "passive tracer"
        elif whatToPlot == "vmag":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "velocity magnitude"
        elif whatToPlot == "p":
            cs = ax.tricontourf(triang, var - p_0, levels = clevels)
            title = "pressure perturbation"
        else:
            sys.exit("Not a valid thing to  plot.")
        fig.colorbar(cs)
        plt.axis('equal')
        plt.axis(radius * 1.025 * np.array([-1,1,-1,1]))

        plt.plot(radius*np.cos(th), radius*np.sin(th), 'k')

        # plt.plot(x[bb], y[bb], 'k.', markersize = 5)

        plt.plot(x[ff], y[ff], 'r.', markersize = 5)

        plt.plot(xNon, yNon, 'ko', markersize = 5)
        
        if plotVectors:
            plt.plot(np.vstack((x-u/2, x+u/2)),
                     np.vstack((y-v/2, y+v/2)), 'k')
        if plotTracers:
            plt.plot(xt, yt, 'k.', markersize = 2)

        plt.title(title + ",  t = {0:4f}".format(tSaves[frame]))
        
        fig.savefig(outstr + '{0:04d}'.format(frame) + '.png',
                    bbox_inches = 'tight')
        
        plt.close()

        frame = frame + 1







