#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import twoBodyTool
import numpy as np
from scipy.integrate import odeint
from twoBodyTool import F, sv2coe

# Defining initial value of position and velocity vector
x0 = -6045
y0 = -3490
z0 = 2500
vx0 = -3.457
vy0 = 6.618
vz0 = 2.533

# Defining period of the orbit in second and its interval
t = np.linspace(0,8200,1000)
# Defining the solution as the result of second ode integral
solution = odeint(F,[x0,y0,z0,vx0,vy0,vz0],t)

# Defining results as aa empty array for plotting
results = []
# Example of loop function of the results data set
for data in range(0,1000):
    r_test = np.squeeze(solution[data:data+1,0:3])
    v_test = np.squeeze(solution[data:data+1,3:6])
    t1 = np.linspace(0,8200,data+1)
    results.append(sv2coe(r_test, v_test))

# Final form of results in array
results = np.array(results)

# Example of how to plot the relationships between orbital elements and the period
import matplotlib.pyplot as plt
x1 = results[:,0]
y1 = t1
plt.plot(y1,x1, '.')
plt.grid(True)
plt.xlabel("Period of orbit (in second)")
plt.ylabel("Semimajor axis value")
plt.title("Variation of Semimajor Axis value for each period of orbital time")
plt.ylim((8785, 8790))
plt.xticks(np.arange(0, 8200, 400))
plt.show()

# In[ ]:




