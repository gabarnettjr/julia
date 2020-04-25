import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

###########################################################################

# Choose the contour levels you want to use
# clevels = np.linspace(-1, 1, 21)
clevels = np.arange( start = -1.05, stop = 1.15, step = .1)

###########################################################################

# Load the x-coordinates of the scattered nodes
x = np.array([])
with open('./results/x.txt') as f:
    for line in f:
        x = np.hstack((x, np.float(line)))

# Load the y-coordinates of the scattered nodes
y = np.array([])
with open('./results/y.txt') as f:
    for line in f:
        y = np.hstack((y, np.float(line)))

# Load the x-coordinates of the regular nodes
xe = np.array([])
with open('./results/xe.txt') as f:
    for line in f:
        xe = np.hstack((xe, np.float(line)))

# Load the y-coordinates of the regular nodes
ye = np.array([])
with open('./results/ye.txt') as f:
    for line in f:
        ye = np.hstack((ye, np.float(line)))

###########################################################################

# plot the nodes and save them
fig = plt.figure(figsize = (12, 10))
plt.plot(xe, ye, '.', x, y, '.')
plt.axis('equal')
plt.title('{0:0d} total nodes'.format(len(x)))
fig.savefig('./figures/nodes.png', bbox_inches = 'tight')

###########################################################################

# Load the approximate solution at the regular nodes
app = np.array([])
with open('./results/app.txt') as f:
    for line in f:
        app = np.hstack((app, np.float(line)))

# Load the exact solution at the regular nodes
exact = np.array([])
with open('./results/exact.txt') as f:
    for line in f:
        exact = np.hstack((exact, np.float(line)))

###########################################################################

# Get the triangular mesh for plotting contours
triang = mtri.Triangulation(xe, ye)

# Contour plot of the approximation
fig = plt.figure(figsize = (12, 10))
ax = fig.add_subplot(111)
cs = ax.tricontourf(triang, app, levels = clevels)
fig.colorbar(cs)
plt.axis('equal')
plt.axis([-1,1,-1,1])
fig.savefig('./figures/app.png', bbox_inches = 'tight')

# Contour plot of the error
fig = plt.figure(figsize = (12, 10))
ax = fig.add_subplot(111)
cs = ax.tricontourf(triang, app-exact, levels = 20)
fig.colorbar(cs)
plt.axis('equal')
plt.axis([-1,1,-1,1])
fig.savefig('./figures/diff.png', bbox_inches = 'tight')

###########################################################################