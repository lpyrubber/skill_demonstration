import numpy as np
import matplotlib.pyplot as plt

m=3
data = np.loadtxt("data.txt")
cp = np.loadtxt("cp.txt")
re = np.loadtxt("residual.txt")
rem = np.loadtxt("residual_map.txt")
#print(data)
x = np.array(data[0,:]).reshape((51,51))
y = np.array(data[1,:]).reshape((51,51))
#print(re[2,:])
O = np.array(data[2,:]).reshape((51,51))
F = np.array(data[4,:]).reshape((51,51))
Rem = np.array(rem[1,:]).reshape((51,51))
R = np.array(re[0,:])
plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#for i in range(51):
#    plt.plot(x[:,i],y[:,i],'-k',linewidth=0.5)
#    plt.plot(x[i,:],y[i,:],'-k',linewidth=0.5)
#levels= np.linspace(np.min(Rem),np.max(Rem),20)
#print(np.min(Rem))
#print(np.max(Rem))
#plt.plot(x[0,:],cp[3,:],linewidth=0.5)
#plt.contour(x,y,F, levels=levels)
#plt.contour(x,y,Rem, levels=levels);
#ax.plot_surface(x,y,Rem,edgecolor='royalblue', lw=0.5, rstride=8, cstride=8,
#                alpha=0.3)
#plt.xlim([-0.1,1.1])
#plt.ylim([-0.3,0.2])
#plt.colorbar()
#plt.title("Mesh near airfoil")
plt.yscale("log")
plt.plot(re[0,:],linewidth=0.5)
plt.plot(re[1,:],linewidth=0.5)
plt.plot(re[2,:],linewidth=0.5)
plt.plot(re[3,:],linewidth=0.5)
plt.plot(re[4,:],linewidth=0.5)
plt.show()
