import numpy as np
import matplotlib.pyplot as plt
import math

data=np.loadtxt("small_cpd.txt",skiprows=1)

l=np.loadtxt("clusters.txt",dtype='int')
lv=np.loadtxt("clusters_v.txt",dtype='int')

m=np.loadtxt("medoids_v.txt")
mv=np.loadtxt("medoids_server.txt")


dim=data.shape[1]
N_points=data.shape[0]
#n_c=l.max()+1
n_c=512
print("N_point=%d, dim=%d, n_c=%d\n"%(N_points,dim,n_c))



for i in range(n_c):
    temp=0
    for j in range(dim):
        temp=temp+(m[i][j]-mv[i][j])*(m[i][j]-mv[i][j])

    print("%d: %e\n"%(i,temp))

"""
#for initial part, the c[i]=i
for i in range(N_points):
    min_v=1e16
    id_v=-1
    for j in range(n_c):
        dis_v=0
        for k in range(dim):
            delta=data[j][k]-data[i][k]
            dis_v=dis_v+delta*delta
        dis_v=math.sqrt(dis_v)
        if(dis_v<min_v):
            min_v=dis_v
            id_v=j
    if(id_v != l[i]):
        print("my part miscoverage at %d (%d<-%d)\n"%(i, l[i],id_v))
    if(id_v !=lv[i]):
        print("TA part miscoverage at %d (%d<-%d)\n"%(i, l[i],id_v))
        

"""