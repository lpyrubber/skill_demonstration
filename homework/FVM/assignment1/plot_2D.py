import numpy as np
import matplotlib.pyplot as plt

method=5;
method_name=['pjc','pgs','ljc','lgs','adl']
titles=['Method 1 Point Jocobi', 'Method 2 Point Gauss-Siedel', 'Method 3 Line Jocobi', 'Method 4 Point Gauss-Siedel', 'Method 5 Alternating Direction Implicit']
data = np.loadtxt("2d_data.txt")
cp = np.loadtxt("2d_cp.txt")
re = np.loadtxt("2d_residual.txt")

x = np.array(data[0,:]).reshape((51,51))
y = np.array(data[1,:]).reshape((51,51))

for i in range(method):
    plt.figure(2*i)
    plt.plot(x[0,:],cp[i,:],'k')
    plt.title(titles[i])
    plt.xlabel("x")
    plt.ylabel("$C_p$")
    plt.xlim([-0.1,1.1])
    plt.ylim([-0.2,0.25])
    plt.savefig(method_name[i]+'_cp.png')

    plt.figure(2*i+1)
    plt.plot(re[i,:],'k')
    plt.title(titles[i])
    plt.xlabel("step")
    plt.ylabel("residual")
    plt.yscale("log")
    plt.savefig(method_name[i]+'_residual.png')

plt.figure(2*method)
plt.yscale("log")
plt.plot(re[0,:],'k',label="pjc")
plt.plot(re[1,:],'--k',label="pgs")
plt.plot(re[2,:],'-.k',label="ljc")
plt.plot(re[3,:],':k',label="lgs")
plt.plot(re[4,:],'0.5',label="adi")
plt.title("Residual For All Methods")
plt.xlabel("step")
plt.ylabel("residual")
plt.legend()
plt.savefig("residual_all.png")

plt.figure(2*method+1)
plt.plot(x[0,:],cp[0,:],'k',label="pjc")
plt.plot(x[0,:],cp[1,:],'--k',label="pgs")
plt.plot(x[0,:],cp[2,:],'0.8',label="ljc",linewidth=5)
plt.plot(x[0,:],cp[3,:],':k',label="lgs",linewidth=4)
plt.plot(x[0,:],cp[4,:],'0.5',label="adi",linewidth=2)
plt.xlim([-0.1,1.1])
plt.ylim([-0.2,0.25])
plt.title("$C_p$ For All Methods")
plt.xlabel("x")
plt.ylabel("$C_p$")
plt.legend()
plt.savefig("cp_all.png")

plt.figure(2*method+2)
for i in range(51):
    plt.plot(x[:,i],y[:,i],'-k',linewidth=0.5)
    plt.plot(x[i,:],y[i,:],'-k',linewidth=0.5)
plt.title("Grid Near Airfoil")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim([-1,2])
plt.ylim([0,3])
plt.savefig("grid.png")
