import numpy as np
import matplotlib.pyplot as plt

CFL=1.0
dt=CFL*0.05
C=1
step=10
time=dt*step
li=0.5
lf=li+time*C
x=np.arange(0,2.05,0.05)
titles = ["Method 1 Explicit Backward", "Method 2 Explicit Forward", "Method 3 Explicit Central", 
          "Method 4 Implicit Central", "Method 5 Crank Nicolson", "Method 6 Lax", "Method 7 Lax Wendroff", 
            "Method 8 MacCormack", "Method 9 Jameson", "Method 10 Warming Beam", "Method 11 Upwind"]
lut_e = [0,1,2,5,6,7,8,9,10]
lut_i = [3,4]

data1 = np.loadtxt("explicit_data.txt")
data2 = np.loadtxt("implicit_data.txt")
af=[0]*41
ai=[0]*41
for i in range(41):
    if (x[i]<=lf):
        af[i]=1.0
    else:
        af[i]=0.5

    if (x[i]<=li):
        ai[i]=1.0
    else:
        ai[i]=0.5
print(af)

print(x)
i=0
for i in range(11):
    if(i<9):
        plt.figure(i)
        
        y=data1[i+2,:]
        plt.plot(x,y,'-k')
        plt.plot(x,ai,'--k',label='Initial condition')
        plt.plot(x,af,'k',linestyle='dashdot', label='Exact solution')
        plt.legend()
#print(titles[lut_e[i]].replace(" ","_"))
        plt.title(titles[lut_e[i]]+' w/ CFL='+str(CFL))
        plt.savefig(titles[lut_e[i]].replace(" ","_")+'_CFL='+str(CFL)+'.png')
    else:
        plt.figure(i)
        i2=i-9
        y=data2[i2+2,:]
        plt.plot(x,y,'-k')
        plt.plot(x,ai,'--k',label='Initial condition')
        plt.plot(x,af,'k',linestyle='dashdot', label='Exact solution')
        plt.legend()
        plt.title(titles[lut_i[i2]]+' w/ CFL='+str(CFL))
        plt.savefig(titles[lut_i[i2]].replace(" ","_")+'_CFL='+str(CFL)+'.png')
