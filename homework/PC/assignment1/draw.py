import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.spatial import Voronoi, voronoi_plot_2d

data=np.loadtxt("tenth_large_cpd.txt",skiprows=1)
#data=np.loadtxt("data.txt",skiprows=1)
#label=np.loadtxt("label.txt")
label=np.loadtxt("clusters.txt",dtype='int')
#medroid=np.loadtxt("medroid.txt");
medroid=np.loadtxt("centroids.txt");


#vor = Voronoi(medroid)
#x=data[:,0]
#y=data[:,1]
index=19
dim=data.shape[1]
N_points=data.shape[0]

n_c=label.max()+1
print("N_point=%d, dim=%d, n_c=%d\n"%(N_points,dim,n_c))
num=np.zeros(n_c)
average=np.zeros((n_c,dim))
error=np.zeros((n_c,2))
for i in range(N_points):
    num[label[i]]=num[label[i]]+1
    temp=0
    for j in range(dim):
        temp=temp+(data[i][j]-medroid[label[i]][j])*(data[i][j]-medroid[label[i]][j])
        average[label[i]][j]=average[label[i]][j]+data[i][j]
    temp=math.sqrt(temp)
    error[label[i]][0]=error[label[i]][0]+temp
for i in range(n_c):
    if (num[i]):
        error[i][0]=error[i][0]/num[i]
    else:
        print("zero at %d, num=%d average=%f\n"%(i,num[i],average[i][0]))
    temp=0
    for j in range(dim):
        average[i][j]=average[i][j]/num[i]
        temp=temp+(average[i][j]-medroid[i][j])*(average[i][j]-medroid[i][j])
    temp=math.sqrt(temp)
    error[i][1]=error[i][1]+temp
#print("medroid location\n")
#print(medroid)
#print("average distance inside group\n")
#print(error[:,0])
#print("average\n")
#print(average)
#print("distance between average and median\n")
#print(error[:,1])
f = open("result_%d.txt"%index,"w") 
#f.write("medroid location\n")
#f.write(" ".join(map(str, medroid)))
f.write("\n\naverage distance inside group\n")
f.write(" ".join(map(str, error[:,0])))

#f.write("\n\naverage\n")
#f.write(" ".join(map(str, average)))

f.write("\n\ndistance between average and median\n")
f.write(" ".join(map(str, error[:,1])))

f.close()
#fig, ax = plt.subplots()
#plt.scatter(x,y)
#fig=voronoi_plot_2d(vor)
#fig=plt.scatter(x,y,c=label)
#fig=plt.scatter(medroid[:,0],medroid[:,1])
#fig=voronoi_plot_2d(vor)
#plt.show()
