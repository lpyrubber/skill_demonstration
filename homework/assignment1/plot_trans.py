import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm


m=0.908
data = np.loadtxt("trans_data.txt")
cp = np.loadtxt("trans_cp.txt")
mach_data = np.loadtxt("trans_mach.txt")
re = np.loadtxt("trans_residual.txt")
x = np.array(data[0,:]).reshape((51,51))
y = np.array(data[1,:]).reshape((51,51))
mach = np.array(mach_data[:,0]).reshape((51,51))
u =    np.array(mach_data[:,1])
v = np.array(mach_data[:,2])
u = u.reshape((51,51))
v = v.reshape((51,51))




plt.figure(0)
plt.plot(x[0,:],cp,'k')
plt.title("Murman-Cole $C_p$ w/ $M_{inf}$="+str(m))
plt.xlabel("x")
plt.ylabel("$C_p$")
plt.xlim([-0.1,1.1])
#for 0.908
plt.ylim([-0.5,0.8])
#for 0.908
#plt.ylim([-0.3,0.3])
plt.savefig('Murman_'+str(m)+'cp.png')

plt.figure(1)
plt.plot(re,'k')
plt.title("Murman-Cole residual w/ $M_{inf}$="+str(m))
plt.xlabel("step")
plt.ylabel("residual")
plt.yscale("log")
plt.savefig('Murman_'+str(m)+'_residual.png')

fig, ax = plt.subplots()
ax.quiver(x,y,u,v,scale=10, width=0.0005, headwidth=50, headlength=50)
ax.set_xlim
ax.set_xlim([-1,2])
#for 0.908
ax.set_ylim([0.0, 3])
#for 0.735
#ax.set_ylim([0.0, 1])
ax.set_title("Mach contours w/ $M_{inf}$="+str(m))
#for 0.908
CS = ax.contour(x, y, mach,levels=[1.0,1.1,1.2,1.3],linewidths=[2.0,2.0,2.0,2.0], colors=['k','k','k','k'])
#for 0.735
#CS = ax.contour(x, y, mach,levels=[0.65,0.7,0.75,0.8],linewidths=[2.0,2.0,2.0,2.0], colors=['k','k','k','k'])
ax.clabel(CS, inline=True, fontsize=12)
plt.savefig('Murman_'+str(m)+'_mach.png')