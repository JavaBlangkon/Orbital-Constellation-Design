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

plt.figure(figsize=(11.69, 16.53 ), dpi=600)

x1 = results[:,0]
y1 = t1
plt.subplot(3, 2, 1)
plt.plot(y1,x1, '.')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel("a (in km)")
plt.title("Variation of Semimajor Axis value for one orbital period")
plt.ylim((8785, 8790))
plt.xticks(np.arange(0, 8200, 1000))

x2 = results[:,1]
y2 = t1
plt.subplot(3, 2, 2)
plt.plot(y2,x2, '.', color='green')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel("i")
plt.title("Variation of Inclination of orbit for one orbital period")
plt.ylim((150, 155))
plt.xticks(np.arange(0, 8200, 1000))

x3 = results[:,2]
y3 = t1
plt.subplot(3, 2, 3)
plt.plot(y3,x3, '.', color='yellow')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel(r'$\Omega$' + ' (in degree)')
plt.title("Variation of Right ascencion of the ascending node for one orbital period")
plt.ylim((255, 256))
plt.xticks(np.arange(0, 8200, 1000))

x4 = results[:,3]
y4 = t1
plt.subplot(3, 2, 4)
plt.plot(y4,x4, '.', color='brown')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel("e")
plt.title("Variation of Eccentricity for one orbital period")
plt.ylim((0, 0.3))
#plt.axis('equal')
plt.xticks(np.arange(0, 8200, 1000))

x5 = results[:,4]
y5 = t1
plt.subplot(3, 2, 5)
plt.plot(y5,x5, '.', color='orange')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel(r'$\omega$' + ' (in degree)')
plt.title("Variation of Argument of perigee for one orbital period")
plt.ylim((10, 25))
plt.xticks(np.arange(0, 8200, 1000))

x6 = results[:,5]
y6 = t1
plt.subplot(3, 2, 6)
plt.plot(y6,x6, '.', color='pink')
plt.grid(True)
plt.xlabel("Time (in second)")
plt.ylabel(r'$\theta$' + ' (in degree)')
plt.title("Variation of True Anomaly for one orbital period")
#plt.ylim((10, 35))
plt.subplots_adjust(hspace=0.3, wspace=0.5)
plt.suptitle("Variations of Classical Orbital Elements over one orbital period", fontsize=18)
#plt.figure(figsize=(11.69, 16.53), dpi=600)
#plt.savefig('coeSimulation.pdf', dpi=600)
plt.savefig('coeSimulationNew.pdf')
plt.show()

# In[ ]:




