import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("./data_1.txt", delimiter = " ")
x=data[:,0]
y=data[:,3]
plt.plot(x,y)
plt.show()
