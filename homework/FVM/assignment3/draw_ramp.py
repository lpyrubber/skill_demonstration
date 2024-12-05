import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm


data = np.loadtxt("cylinder_points.txt")
#data = np.loadtxt("ramp_points.txt")
slist = np.loadtxt("slist.txt")
clist = np.loadtxt("clist.txt")
plist = np.loadtxt("plist.txt")
cell = np.loadtxt("cell_info.txt")

Np=data.shape[0]
Ns=slist.shape[0]
Nc=clist.shape[0]

print("Np=%d, Ns=%d, Nc=%d\n" % (Np, Ns, Nc))

x=np.zeros(2)
y=np.zeros(2)

#print("slist=%d %d, xp = %f %f\n" % (slist[1][0],slist[1][1],data[0][0],data[0][1]))

plt.figure(0)
for i in range(Ns):
    p1=int(slist[i][0])
    p2=int(slist[i][1])
    x[0]=data[p1][0]
    x[1]=data[p2][0]
    y[0]=data[p1][1]
    y[1]=data[p2][1]
    plt.plot(x,y,'k')
#ax.contour(xc,yc,zc)
#plt.scatter(data[:,0],data[:,1])
plt.scatter(cell[:,1],cell[:,2])
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(False)  # Disable default grid to match the style
plt.title("Cylinder Grid 80x80 with cell center")
plt.xlabel('x')
plt.ylabel('y')
plt.show()