import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


r_rho=pd.read_csv('rho_exact.csv',sep=",",header=None)
r_p=pd.read_csv('p_exact.csv',sep=",",header=None)
r_u=pd.read_csv('u_exact.csv',sep=",",header=None)
data=np.loadtxt("data.txt")
x=data[:,0]
y=data[:,1]
z=data[:,2]
w=data[:,3]
plt.figure(0)
plt.title('Density distribution at t=0.5')
plt.plot(r_rho[0],r_rho[1],'k',linestyle='dashdot', label='Exact solution')
plt.plot(x,y,'k')
plt.legend()
plt.savefig('density.png')

plt.figure(1)
plt.title('Velocity distribution at t=0.5')
plt.plot(r_u[0],r_u[1],'k',linestyle='dashdot', label='Exact solution')
plt.plot(x,z,'k')
plt.legend()
plt.savefig('velocity.png')


plt.figure(2)
plt.title('Pressure distribution at t=0.5')
plt.plot(r_p[0],r_p[1],'k',linestyle='dashdot', label='Exact solution')
plt.plot(x,w,'k')
plt.legend()
plt.savefig('pressure.png')

