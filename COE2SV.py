#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

#G_const for gravitational constant, ME for mass of earth, RE for radius of earth
G_const = 6.67408 * 10**(-11)
ME = 5.972364730419773 * 10**24
RE = 6378100
#Declaring value of miu
miu = (G_const * ME) * 10**(-9)

# defining function coe2sv(h, i, omega_capt, e, omega_case, theta) for transforming state vector to classical orbital elements
# in the following definition of function sv2coe, sv is a state vector and coe is classical orbital elements
# the input element of the function sv2coe that need to be defined consists of:
#   h = the magnitude of specific angular momentum (in km^2/s)
#   i = the magnitude of inclination (in degree)
#   omega_capt = the magnitude of right ascencion of the ascending node (in degree)
#   e = the magnitude of eccentricity (no unit)
#   omega_case = the magnitude of argument of perigee (in degree)
#   theta = the magnitude of true anomaly (in degree)
#   all the inputs should be in the correct unit
# 
# the possible output of this function are as follows:
#   r_peri = the position vector in perifocal coordinates (in km)
#   v_peri = the velocity vector in perifocal coordinates (in km/s)
#   Q_geo = the matrix Q_geo transformation from perifocal to geocentric equatorial coordinates (no unit)
#   r_geo = the position vector in geocentric frame (in km)
#   v_geo = the velocity vector in geocentric frame (in km/s)

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

print("The geocentric equatorial position vector is", r_geo, "km")
print("The geocentric equatorial velocity vector is", v_geo, "km/s")


# In[ ]:




