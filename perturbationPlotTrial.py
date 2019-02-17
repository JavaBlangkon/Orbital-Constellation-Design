#!/usr/bin/env python
# coding: utf-8

# In[13]:


import twoBodyTool
import numpy as np
from scipy.integrate import odeint
from twoBodyTool import F, F1, sv2coe

x0 = -6045
y0 = -3490
z0 = 2500
vx0 = -3.457
vy0 = 6.618
vz0 = 2.533

t = np.linspace(0,8200 * 11, 100000)

solution1 = odeint(F1, [x0, y0, z0, vx0, vy0, vz0], t)

results1 = []

for data in range(0, 100000):
    r_test = np.squeeze(solution1[data:data+1, 0:3])
    v_test = np.squeeze(solution1[data:data+1, 3:6])
    t1 = np.linspace(0,(8200 * 11) / 3600, data+1)
    results1.append(sv2coe(r_test, v_test))
    
results1 = np.array(results1)

# a = e = i = 0
# [a, i, omega_capt, e, omega_case, theta, h, T]


# In[8]:


import matplotlib.pyplot as plt
x = solution1[:, 0]
y = solution1[:, 1]
plt.plot(x, y, '.')
plt.axis('equal')


# In[19]:


import matplotlib.pyplot as plt

plt.figure(figsize=(16.53, 11.69), dpi=600)

x1 = results1[:, 0]
y1 = t1
plt.subplot(3, 2, 1)
plt.plot(y1, x1, '.')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\alpha$' + ' (in km)')
plt.title("Variation of Semimajor Axis for one orbital period \n")
#plt.ylim((8785, 8790))
#plt.xticks(np.arange(0, 8200, 1000))

x2 = results1[:, 1]
y2 = t1
plt.subplot(3, 2, 2)
plt.plot(y2, x2, '.', color='green')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel('i' + ' (in degree)')
plt.title("Variation of Inclination of orbit for one orbital period \n")
#plt.ylim((150, 155))
#plt.xticks(np.arange(0, 8200, 1000))

x3 = results1[:, 2]
y3 = t1
plt.subplot(3, 2, 3)
plt.plot(y3, x3, '.', color='yellow')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\Omega$' + ' (in degree)')
plt.title("Variation of Right ascencion of the ascending node for one orbital period \n")
#plt.ylim((255, 256))
#plt.xticks(np.arange(0, 8200, 1000))

x4 = results1[:, 3]
y4 = t1
plt.subplot(3, 2, 4)
plt.plot(y4, x4, '.', color='brown')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel("e")
plt.title("Variation of Eccentricity for one orbital period \n")
#plt.ylim((0, 0.3))
#plt.axis('equal')
#plt.xticks(np.arange(0, 8200, 1000))

x5 = results1[:, 4]
y5 = t1
plt.subplot(3, 2, 5)
plt.plot(y5, x5, '.', color='orange')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\omega$' + ' (in degree)')
plt.title("Variation of Argument of perigee for one orbital period \n")
#plt.ylim((10, 25))
#plt.xticks(np.arange(0, 8200, 1000))

x6 = results1[:, 5]
y6 = t1
plt.subplot(3, 2, 6)
plt.plot(y6, x6, '.', color='pink')
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\theta$' + ' (in degree)')
plt.title("Variation of True Anomaly for one orbital period \n")
#plt.ylim((10, 35))
plt.subplots_adjust(hspace=0.5, wspace=0.3)
plt.suptitle("Variations of Classical Orbital Elements over one orbital period \n under J2 perturbation", fontsize=18, fontweight=700)
#plt.figure(figsize=(11.69, 16.53), dpi=600)
#plt.savefig('coeSimulation.pdf', dpi=600)
plt.savefig('coeSimulationNewPerturbation.pdf')
plt.show()


# In[ ]:




