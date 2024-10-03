import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt("data.txt",skiprows=1)
label=np.loadtxt("label.txt")
medroid=np.loadtxt("medroid.txt");
print(label)

x=data[:,0]
y=data[:,1]

#plt.scatter(x,y)
plt.scatter(x,y,c=label)
plt.scatter(medroid[:,0],medroid[:,1])
plt.show()
