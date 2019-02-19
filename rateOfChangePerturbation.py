#!/usr/bin/env python
# coding: utf-8

# In[1]:


import twoBodyTool
import numpy as np
from scipy.integrate import odeint
from twoBodyTool import F, F1, sv2coe, RE, J2, miu

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

rateFactor = results1[:1,8]*J2*((RE*10**(-3)/results1[:1,9])**2)

omgRate = -(3/2)*(rateFactor*np.cos(results1[:1,1]*np.pi/180))
omgCaseRate = (3/4)*rateFactor*(5*(np.cos(results1[:1,1]*np.pi/180)**2)-1)


# In[9]:


import numpy as np
import pylab as plot
import matplotlib.pyplot as plt
import numpy, scipy, pylab, random
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from scipy import stats


x1 = results1[:, 4]
y1 = t1

slope, intercept, r_value, p_value, std_err = stats.linregress(y1,x1)
line = slope*y1+intercept
plt.plot(y1, line, 'red', label='fitted line')

plt.scatter(y1,x1,color='green', s=3)

plt.show()


# In[10]:


import numpy as np
import pylab as plot
import matplotlib.pyplot as plt
import numpy, scipy, pylab, random
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from scipy import stats

plt.figure(figsize=(16.53, 11.69), dpi=600)

x1 = t1
y1 = results1[:, 0]
plt.subplot(3, 2, 1)
slope, intercept, r_value, p_value, std_err = stats.linregress(x1,y1)
line = slope*x1+intercept
plt.plot(x1, line, 'red', label='fitted line')
plt.scatter(x1,y1,color='blue', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\alpha$' + ' (in km)')
plt.title("Variation of Semimajor Axis for one orbital period \n")
#plt.ylim((8785, 8790))
#plt.xticks(np.arange(0, 8200, 1000))

x2 = t1
y2 = results1[:, 1]
plt.subplot(3, 2, 2)
slope, intercept, r_value, p_value, std_err = stats.linregress(x2,y2)
line = slope*x2+intercept
plt.plot(x2, line, 'red', label='fitted line')
plt.scatter(x2,y2,color='green', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel('i' + ' (in degree)')
plt.title("Variation of Inclination of orbit for one orbital period \n")
#plt.ylim((150, 155))
#plt.xticks(np.arange(0, 8200, 1000))

x3 = t1
y3 = results1[:, 2]
plt.subplot(3, 2, 3)
slope, intercept, r_value, p_value, std_err = stats.linregress(x3,y3)
line = slope*x3+intercept
plt.plot(x3, line, 'red', label='fitted line')
plt.scatter(x3,y3,color='yellow', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\Omega$' + ' (in degree)')
plt.title("Variation of Right ascencion of the ascending node for one orbital period \n")
#plt.ylim((255, 256))
#plt.xticks(np.arange(0, 8200, 1000))

x4 = t1
y4 = results1[:, 3]
plt.subplot(3, 2, 4)
slope, intercept, r_value, p_value, std_err = stats.linregress(x4,y4)
line = slope*x4+intercept
plt.plot(x4, line, 'red', label='fitted line')
plt.scatter(x4,y4,color='brown', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel("e")
plt.title("Variation of Eccentricity for one orbital period \n")
#plt.ylim((0, 0.3))
#plt.axis('equal')
#plt.xticks(np.arange(0, 8200, 1000))

x5 = t1
y5 = results1[:, 4]
plt.subplot(3, 2, 5)
slope, intercept, r_value, p_value, std_err = stats.linregress(x5,y5)
line = slope*x5+intercept
plt.plot(x5, line, 'red', label='fitted line')
plt.scatter(x5,y5,color='orange', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\omega$' + ' (in degree)')
plt.title("Variation of Argument of perigee for one orbital period \n")
#plt.ylim((10, 25))
#plt.xticks(np.arange(0, 8200, 1000))

x6 = t1
y6 = results1[:, 5]
plt.subplot(3, 2, 6)
slope, intercept, r_value, p_value, std_err = stats.linregress(x6,y6)
line = slope*x6+intercept
plt.plot(x6, line, 'red', label='fitted line')
plt.scatter(x6,y6,color='pink', s=3)
plt.grid(True)
plt.xlabel("Time (in hours)")
plt.ylabel(r'$\theta$' + ' (in degree)')
plt.title("Variation of True Anomaly for one orbital period \n")
#plt.ylim((10, 35))
plt.subplots_adjust(hspace=0.5, wspace=0.3)
plt.suptitle("Variations of Classical Orbital Elements over one orbital period \n under J2 perturbation", fontsize=18, fontweight=700)
#plt.figure(figsize=(11.69, 16.53), dpi=600)
#plt.savefig('coeSimulation.pdf', dpi=600)
plt.savefig('regressionPerturbation.pdf')
plt.show()


# In[ ]:




