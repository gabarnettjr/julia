 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os

#####################################################################

# USER INPUT

autoLevels = True
# clevels = np.arange(0, 21, 1)
# clevels = np.arange(-.025, 1.075, .05)

# It can automatically choose the contour levels if this is True
if autoLevels:
    clevels = 20

# If the user passes in 3 arguments, assume the first argument is whether
# or not to plot vectors (0 or 1), the second argument is whether or not to
# plot moving points (0 or 1), and the third argument is what to plot, which
# can be vmag, u, v, rho, q, or p.
if len(sys.argv) >= 4:
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
    plotTracers = True
    whatToPlot = "vmag"

# Directory where the *.txt files are located (in string)
if len(sys.argv) == 5:
    instr = sys.argv[4]
else:
    instr = 'results/'

# Get the variable string
varstr = instr + whatToPlot

# Directory where the *.png files will end up (out string)
outstr = 'figures/'

#####################################################################

# Load the constants and background states

Cp = np.array([], float)
with open(instr + 'Cp.txt') as f:
    for line in f:
        Cp = np.hstack((Cp, float(line)))
Cp = Cp[0]

Cv = np.array([], float)
with open(instr + 'Cv.txt') as f:
    for line in f:
        Cv = np.hstack((Cv, float(line)))
Cv = Cv[0]

R = np.array([], float)
with open(instr + 'R.txt') as f:
    for line in f:
        R = np.hstack((R, float(line)))
R = R[0]

pi_0 = np.array([], float)
with open(instr + 'pi_0.txt') as f:
    for line in f:
        pi_0 = np.hstack((pi_0, float(line)))

th_0 = np.array([], float)
with open(instr + 'th_0.txt') as f:
    for line in f:
        th_0 = np.hstack((th_0, float(line)))

T_0 = np.array([], float)
with open(instr + 'T_0.txt') as f:
    for line in f:
        T_0 = np.hstack((T_0, float(line)))

###############################################################################

if autoLevels:

    minPi = np.array([], float)
    with open(instr + 'minPi.txt') as f:
        for line in f:
            minPi = np.hstack((minPi, float(line)))
    minPi = minPi[0]
    
    maxPi = np.array([], float)
    with open(instr + 'maxPi.txt') as f:
        for line in f:
            maxPi = np.hstack((maxPi, float(line)))
    maxPi = maxPi[0]
    
    minU = np.array([], float)
    with open(instr + 'minU.txt') as f:
        for line in f:
            minU = np.hstack((minU, float(line)))
    minU = minU[0]
    
    maxU = np.array([], float)
    with open(instr + 'maxU.txt') as f:
        for line in f:
            maxU = np.hstack((maxU, float(line)))
    maxU = maxU[0]
    
    minV = np.array([], float)
    with open(instr + 'minV.txt') as f:
        for line in f:
            minV = np.hstack((minV, float(line)))
    minV = minV[0]
    
    maxV = np.array([], float)
    with open(instr + 'maxV.txt') as f:
        for line in f:
            maxV = np.hstack((maxV, float(line)))
    maxV = maxV[0]
    
    minTh = np.array([], float)
    with open(instr + 'minTh.txt') as f:
        for line in f:
            minTh = np.hstack((minTh, float(line)))
    minTh = minTh[0]
    
    maxTh = np.array([], float)
    with open(instr + 'maxTh.txt') as f:
        for line in f:
            maxTh = np.hstack((maxTh, float(line)))
    maxTh = maxTh[0]

    minQ = np.array([], float)
    with open(instr + 'minQ.txt') as f:
        for line in f:
            minQ = np.hstack((minQ, float(line)))
    minQ = minQ[0]
    
    maxQ = np.array([], float)
    with open(instr + 'maxQ.txt') as f:
        for line in f:
            maxQ = np.hstack((maxQ, float(line)))
    maxQ = maxQ[0]
    
    minT = np.array([], float)
    with open(instr + 'minT.txt') as f:
        for line in f:
            minT = np.hstack((minT, float(line)))
    minT = minT[0]
    
    maxT = np.array([], float)
    with open(instr + 'maxT.txt') as f:
        for line in f:
            maxT = np.hstack((maxT, float(line)))
    maxT = maxT[0]

    # Load the maximum velocity magnitude over the whole time interval
    maxVel = np.array([], float)
    with open(instr + 'maxVel.txt') as f:
        for line in f:
            maxVel = np.hstack((maxVel, float(line)))
    maxVel = maxVel[0]

    # Set the contour levels
    if whatToPlot == "pi":
        clevels = np.linspace(minPi, maxPi, clevels)
    elif whatToPlot == "u":
        clevels = np.linspace(minU, maxU, clevels)
    elif whatToPlot == "v":
        clevels = np.linspace(minV, maxV, clevels)
    elif whatToPlot == "th":
        clevels = np.linspace(minTh, maxTh, clevels)
    elif whatToPlot == "q":
        clevels = np.linspace(minQ, maxQ, clevels)
    elif whatToPlot == "T":
        clevels = np.linspace(minT, maxT, clevels)
    elif whatToPlot == "vmag":
        clevels = np.linspace(0, maxVel, clevels)
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
        x = np.hstack((x, float(line)))

# Load the y-coordinates of the points
y = np.array([], float)
with open(instr + 'y.txt') as f:
    for line in f:
        y = np.hstack((y, float(line)))

# Load the x-coordinates of the NON points
xNon = np.array([], float)
with open(instr + 'xNon.txt') as f:
    for line in f:
        xNon = np.hstack((xNon, float(line)))

# Load the y-coordinates of the NON points
yNon = np.array([], float)
with open(instr + 'yNon.txt') as f:
    for line in f:
        yNon = np.hstack((yNon, float(line)))

# Load the times that will be plotted:
tSaves = np.array([], float)
with open(instr + 'tSaves.txt') as f:
    for line in f:
        tSaves = np.hstack((tSaves, float(line)))

# Load the approximate node spacing
dr = np.array([], float)
with open(instr + 'dr.txt') as f:
    for line in f:
        dr = np.hstack((dr, float(line)))
dr = dr[0]

###############################################################################

# plot the nodes and save the plot

# index of the boundary nodes
bb = np.array([], int)
with open(instr + 'bb.txt') as f:
    for line in f:
        bb = np.hstack((bb, int(line) - 1))

# index of the "fan" nodes
ff = np.array([], int)
with open(instr + 'ff.txt') as f:
    for line in f:
        ff = np.hstack((ff, int(line) - 1))

fig = plt.figure(figsize = (12, 10))
plt.plot(x, y, 'k.', \
         x[bb], y[bb], 'ks', \
         x[ff], y[ff], 'rs', \
         xNon, yNon, 'kx', \
         markersize = 1, fillstyle = 'none')
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
if whatToPlot == "T":
    pi = np.zeros(np.shape(x))
    th = np.zeros(np.shape(x))
if plotTracers:
    xt = np.zeros(np.shape(x))
    yt = np.zeros(np.shape(x))

# Initialize the frame number
frame = 0

# The main loop that plots and saves figures

while True:

    try:

        if plotVectors or (whatToPlot == "vmag"):
            i = 0
            with open(instr + "u" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    u[i] = float(line)
                    i = i + 1
            i = 0
            with open(instr + "v" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    v[i] = float(line)
                    i = i + 1
            if (whatToPlot == "vmag"):
                var = np.sqrt(u**2 + v**2)
            if plotVectors:
                u = u / maxVel * dr
                v = v / maxVel * dr
        elif whatToPlot == "T":
            i = 0
            with open(instr + "pi" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    pi[i] = float(line)
                    i = i + 1
            i = 0
            with open(instr + "th" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    th[i] = float(line)
                    i = i + 1
            var = pi * th
        else:
            i = 0
            with open(varstr + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    var[i] = float(line)
                    i = i + 1
        if plotTracers:
            i = 0
            with open(instr + "xt" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    xt[i] = float(line)
                    i = i + 1
            i = 0
            with open(instr + "yt" + "_{0:04d}.txt".format(frame)) as f:
                for line in f:
                    yt[i] = float(line)
                    i = i + 1

    except:

        break

    else:
        
        fig = plt.figure(figsize = (12, 10))
        ax = fig.add_subplot(111)
        if whatToPlot == "pi":
            cs = ax.tricontourf(triang, var-pi_0, levels = clevels)
            title = "Exner pressure perturbation"
        elif whatToPlot == "u":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "velocity west-to-east"
        elif whatToPlot == "v":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "velocity south-to-north"
        elif whatToPlot == "th":
            cs = ax.tricontourf(triang, var-th_0, levels = clevels)
            title = "potential temperature perturbation"
        elif whatToPlot == "q":
            cs = ax.tricontourf(triang, var, levels = clevels)
            title = "passive tracer"
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

        plt.plot(x[ff], y[ff], 'r.', markersize = 5)

        plt.plot(xNon, yNon, 'ks', markersize = 4)
        
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







