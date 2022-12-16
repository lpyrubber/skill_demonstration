import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("./data.txt", delimiter = " ")
x=data[:,0]
y=data[:,1]
plt.plot(x,y)
plt.show()
