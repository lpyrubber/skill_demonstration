import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.txt")
#print(data)
x = np.array(data[0,:]).reshape((51,51))
y = np.array(data[1,:]).reshape((51,51))
y = y.transpose()
print(y)

plt.figure()
for i in range(51):
    plt.plot(x[:,i],y[:,i],'-k',linewidth=0.5)
    plt.plot(x[i,:],y[i,:],'-k',linewidth=0.5)

#plt.xlim([-1,2])
#plt.ylim([0,1.5])
plt.title("Mesh near airfoil")
plt.savefig('grid.png')
