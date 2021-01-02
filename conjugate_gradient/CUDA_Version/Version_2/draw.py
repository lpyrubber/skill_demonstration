import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("./data.txt", delimiter = " ")
plt.plot(data)
plt.show()
