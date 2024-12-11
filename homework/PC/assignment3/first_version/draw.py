import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.spatial import Voronoi, voronoi_plot_2d

#data=np.loadtxt("small_gaussian.txt",skiprows=1)
data=np.loadtxt("data.txt",skiprows=1)
#label=np.loadtxt("label.txt")
label=np.loadtxt("clusters.txt")
#medroid=np.loadtxt("medroid.txt");
medroid=np.loadtxt("centroids.txt");


#vor = Voronoi(medroid)
x=data[:,0]
y=data[:,1]
n_c=int(label.max()+1)
num=np.zeros(n_c)
average=np.zeros((n_c,2))
error=np.zeros((n_c,2))
for i in range(len(data)):
    num[int(label[i])]=num[int(label[i])]+1
    error[int(label[i])][0]=error[int(label[i])][0]+math.sqrt((data[i][0]-medroid[int(label[i])][0])*(data[i][0]-medroid[int(label[i])][0])+(data[i][1]-medroid[int(label[i])][1])*(data[i][1]-medroid[int(label[i])][1]))
    average[int(label[i])][0]=average[int(label[i])][0]+data[i][0]
    average[int(label[i])][1]=average[int(label[i])][1]+data[i][1]
for i in range(n_c):
    error[i][0]=error[i][0]/num[i]
    average[i][0]=average[i][0]/num[i]
    average[i][1]=average[i][1]/num[i]
    error[i][1]=math.sqrt((average[i][0]-medroid[i][0])*(average[i][0]-medroid[i][0])+(average[i][1]-medroid[i][1])*(average[i][1]-medroid[i][1]))
print("medroid location\n")
print(medroid)
print("average distance inside group\n")
print(error[:,0])
print("average\n")
print(average)
print("distance between average and median\n")
print(error[:,1])
#fig, ax = plt.subplots()
#plt.scatter(x,y)
#fig=voronoi_plot_2d(vor)
fig=plt.scatter(x,y,c=label)
fig=plt.scatter(medroid[:,0],medroid[:,1])
#fig=voronoi_plot_2d(vor)
plt.show()
