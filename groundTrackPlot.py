#!/usr/bin/env python
# coding: utf-8

# In[5]:


import twoBodyTool
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from twoBodyTool import F, F1, sv2coe, RE, J2, miu, rate, sphericalCoor

# state vector of satellite (in this case ISS)
x0 = 5720.95209152
y0 = -3642.08487137
z0 = 10.20957305
vx0 = 2.54734135
vy0 = 4.01805211
vz0 = 6.01104338

# period of satellite revolution, with 10000 samples generated
t = np.linspace(0,5544, 100)

# finding the gravitational force under J2 perturbation
solution1 = odeint(F1, [x0, y0, z0, vx0, vy0, vz0], t)

results1 = []
results2 = []

for data in range(0, 100):
    r_test = np.squeeze(solution1[data:data+1, 0:3])
    v_test = np.squeeze(solution1[data:data+1, 3:6])
    t1 = np.linspace(0,(5544)/3600,data+1)
    results1.append(sv2coe(r_test, v_test))
    results2.append(sphericalCoor(r_test))

results1 = np.array(results1)
results2 = np.array(results2)


# In[6]:


from mpl_toolkits.basemap import Basemap

# set up the basic map size and configuration
f = plt.figure(figsize=(10,7.5))
m = Basemap(projection="mill", lon_0=0)

m.drawcoastlines()
m.drawparallels(np.arange(-90,91,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180,181,30), labels=[0,0,0,1])

# defining the longitude and latitude of satellite in array - need more explanation on defining these values
lon = results2[:,1]
lat = results2[:,2]

# assigning longitude and latitude in x and y axis of the map
x,y = m(lon, lat)
# plotting the ground track on 2D map
m.plot(x, y, color="red", latlon=False, marker='.', linestyle='None')


# In[4]:


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

X = solution1[:,0]
Y = solution1[:,1]
Z = solution1[:,2]

scat = ax.scatter(X, Y, Z)

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
    ax.plot([xb], [yb], [zb], 'w')

plt.grid()
plt.show()

