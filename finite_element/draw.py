import numpy as np
import matplotlib.pyplot as plt
a=np.loadtxt('triangle_x.txt')
b=np.loadtxt('triangle_y.txt')
c=np.loadtxt('data_x.txt')
d=np.loadtxt('data_y.txt')
for i in range(len(a)):
    plt.plot([a[i][0],a[i][1]],[b[i][0],b[i][1]],'b')
    plt.plot([a[i][1],a[i][2]],[b[i][1],b[i][2]],'b')
    plt.plot([a[i][2],a[i][0]],[b[i][2],b[i][0]],'b')
    plt.plot([c[i][0],c[i][1]],[d[i][0],d[i][1]],'r')
    plt.plot([c[i][1],c[i][2]],[d[i][1],d[i][2]],'r')
    plt.plot([c[i][2],c[i][0]],[d[i][2],d[i][0]],'r')
    
plt.show()
