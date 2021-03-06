import numpy as np
import matplotlib.pyplot as plt
a=np.loadtxt('triangle_x.txt')
b=np.loadtxt('triangle_y.txt')
c=np.loadtxt('point.txt')
plt.figure()
for i in range(len(a)):
    plt.plot([a[i][0],a[i][1]],[b[i][0],b[i][1]],'k')
    plt.plot([a[i][1],a[i][2]],[b[i][1],b[i][2]],'k')
    plt.plot([a[i][2],a[i][0]],[b[i][2],b[i][0]],'k')
plt.show()

plt.figure()
for i in range(len(c)):
    plt.scatter(c[i][0],c[i][1], c='k', marker = 'o')
plt.show()
