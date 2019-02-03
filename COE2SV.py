#!/usr/bin/env python
# coding: utf-8

# In[1]:


#astropy is a python library for astronomy
#G for gravitational constant, M_earth for mass of earth, R_earth for radius of earth
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import numpy as np

#declaring value of miu and radius of earth
miu = (G.value * M_earth.value) * 10**(-9)
Re = R_earth.value

#defining function coe2sv for transforming classical orbital elements to state vector
def coe2sv(h, i, omega_capt, e, omega_case, theta):
    #Calculating position vector in perifocal coordinates
    r_peri = (h**2 / miu) * (1 / (1 + e * np.cos(theta))) * np.array([[np.cos(theta)], [np.sin(theta)], [0]])
    #Calculating velocity vector in perifocal coordinates
    v_peri = (miu / h) * np.array([[-np.sin(theta)], [e+np.cos(theta)], [0]])
    #Calculating matrix Q_geo transformation from perifocal to geocentric equatorial coordinates
    A = (np.cos(omega_capt) * np.cos(omega_case) - np.sin(omega_capt) * np.sin(omega_case) * np.cos(i))
    B = (-np.cos(omega_capt) * np.sin(omega_case) - np.sin(omega_capt) * np.cos(i) * np.cos(omega_case))
    C = (np.sin(omega_capt) * np.sin(i))
    D = (np.sin(omega_capt) * np.cos(omega_case) + np.cos(omega_capt) * np.cos(i) * np.sin(omega_case))
    E = (-np.sin(omega_capt) * np.sin(omega_case) + np.cos(omega_capt) * np.cos(i) * np.cos(omega_case))
    F = (np.cos(omega_capt) * np.sin(i))
    G = (np.sin(i) * np.sin(omega_case))
    H = (np.sin(i) * np.cos(omega_case))
    I = (np.cos(i))
    Q_geo = np.array([[A, B, C], [D, E, F], [G, H, I]])
    #Transforming perifocal position and velocity to geocentric frame
    r_geo = np.dot(Q_geo, r_peri)
    v_geo = np.dot(Q_geo, v_peri)
    
    return Q_geo, r_peri, v_peri, r_geo, v_geo

#testing using the parameters from Orbital Mechanics page 175 Example 4.5
h = 80000
i = 30 * np.pi / 180
omega_capt = 40 * np.pi / 180
e = 1.4
omega_case = 60 * np.pi / 180
theta = 30 * np.pi / 180

Q_geo, r_peri, v_peri, r_geo, v_geo = coe2sv(h, i, omega_capt, e, omega_case, theta)

print("The geocentric equatorial position vector is", r_geo, u.km)
print("The geocentric equatorial velocity vector is", v_geo, u.km/u.s)


# In[ ]:




