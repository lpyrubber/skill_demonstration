import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

x=np.linspace(-10,2,800)
#print(x)
x_r=pow(10,x)
print(x_r)
y=8.87424*1e-7/x_r
plt.figure(0)

plt.plot(x_r,y)
plt.title("percentage of valid collision in log scale")
plt.xlabel("Kf, reaction rate coefficient unit:mol L-1 s-1 or mol dm-3 s")
plt.ylabel("percentage of valid collision")
plt.yscale('log')
plt.xscale('log')
plt.show()