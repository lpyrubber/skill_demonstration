import numpy as np
import matplotlib.pyplot as plt

m=2
data = np.loadtxt("data.txt")
cp = np.loadtxt("cp.txt")
re = np.loadtxt("residual.txt")
#print(data)
x = np.array(data[0,:]).reshape((51,51))
y = np.array(data[1,:]).reshape((51,51))
print(re[1,:])

O = np.array(data[2,:]).reshape((51,51))
F = np.array(data[3,:]).reshape((51,51))

plt.figure()
#for i in range(51):
#    plt.plot(x[:,i],y[:,i],'-k',linewidth=0.5)
#    plt.plot(x[i,:],y[i,:],'-k',linewidth=0.5)
levels= np.linspace(np.min(F),np.max(F),20)
#plt.plot(x[0,:],cp,linewidth=0.5)
#plt.contour(x,y,F, levels=levels)
#plt.xlim([-0.1,1.1])
#plt.ylim([-0.3,0.2])
plt.title("Mesh near airfoil")

plt.plot(re[0,:],linewidth=0.5)
plt.plot(re[1,:],linewidth=0.5)
plt.savefig('grid.png')
