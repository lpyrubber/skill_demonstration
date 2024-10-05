import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

#data=np.loadtxt("small_gaussian.txt",skiprows=1)
data=np.loadtxt("data.txt",skiprows=1)
#label=np.loadtxt("label.txt")
label=np.loadtxt("clusters.txt")
medroid=np.loadtxt("medroid.txt");
print(medroid)
#vor = Voronoi(medroid)
x=data[:,0]
y=data[:,1]

#fig, ax = plt.subplots()
#plt.scatter(x,y)
#fig=voronoi_plot_2d(vor)
fig=plt.scatter(x,y,c=label)
fig=plt.scatter(medroid[:,0],medroid[:,1])
#fig=voronoi_plot_2d(vor)
plt.show()
