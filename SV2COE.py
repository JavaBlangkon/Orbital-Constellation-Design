#!/usr/bin/env python
# coding: utf-8

# In[23]:


#astropy is a python library for astronomy
#G for gravitational constant, M_earth for mass of earth, R_earth for radius of earth
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import numpy as np

#declaring value of miu and radius of earth
miu = (G.value * M_earth.value) * 10**(-9)
Re = R_earth.value

#declaring unit vector
i_unit = np.array([1, 0, 0])
j_unit = np.array([0, 1, 0])
k_unit = np.array([0, 0, 1])

# defining function sv2coe(r_vec, v_vec) for transforming state vector to classical orbital elements
# in the following definition of function sv2coe, sv is a state vector and coe is classical orbital elements
# the input element of the function sv2coe that need to be defined consists of:
#   r_vec is position vector of the satellite
#   v_vec is velocity vector of the satellite
#   all the input vector should be in 3x1 matrix
# 
# the possible output of this function are as follows:
#   r = the magnitude of position vector (in km)
#   v = the magnitude of velocity vector (in km/s)
#   v_rad = the radial velocity (in km/s)
#   h_vec = the specific angular momentum vector (in km^2/s)
#   h = the magnitude of specific angular momentum (in km^2/s)
#   i = the magnitude of inclination (in degree)
#   N_vec = the node line vector (no unit)
#   N = the magnitude of node line (no unit)
#   omega_capt = the magnitude of right ascencion of the ascending node (in degree)
#   e_vec = the eccentricity vector (no unit)
#   e = the magnitude of eccentricity (no unit)
#   omega_case = the magnitude of argument of perigee (in degree)
#   theta = the magnitude of true anomaly (in degree)
#   rp = the magnitude of radius of perigee (in km)
#   ra = the magnitude of radius of apogee (in km)
#   a = the magnitude of semimajor axis (in km)
#   n = the magnitude of mean motion (no unit)
#   T = the magnitude of period the orbit (in second)

def sv2coe(r_vec, v_vec):
    #Calculating distance from r_vec
    r = np.sqrt(np.dot(r_vec, r_vec))
    #Calculating speed from v_vec
    v = np.sqrt(np.dot(v_vec, v_vec))
    #Calculating radial velocity
    v_rad = (np.dot(r_vec, v_vec)) / r
    if v_rad > 0:
        print("Satellite is flying away from perigee")
    else:
        print("Satellite is flying towards perigee")
    #Calculating specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    #Calculating magnitude of specific angular momentum
    h = np.sqrt(np.dot(h_vec, h_vec)) #first orbital element
    #Calculating inclination
    i = np.arccos(h_vec[2] / h) 
    i = i*180/np.pi #second orbital element
    if i <= 90:
        print("The orbit is retrograde")
    elif 90 < i < 180:
        print("The orbit is posigrade")
    else:
        print("The orbit is greater than 180 and in quadrant ambiguity")
    #Calculating vector node line
    N_vec = np.cross(k_unit, h_vec)
    #Calculating magnitude of vector node line
    N = np.sqrt(np.dot(N_vec, N_vec))
    #Calculating right ascencion of the ascending node
    if N_vec[1] >= 0:
        omega_capt = np.arccos(N_vec[0] / N)
    else:
        omega_capt = 2 * np.pi - np.arccos(N_vec[0] / N) 
    omega_capt = omega_capt * 180 / np.pi #third orbital element
    #Calculating eccentricity vector
    e_vec = (1 / miu) * (np.dot((v**2 - (miu / r)), r_vec) - (np.dot((r * v_rad), v_vec)))
    #Calculating magnitude of eccentricity
    e = np.sqrt(np.dot(e_vec, e_vec)) #fourth orbital element
    #Calculating argument of perigee
    if e_vec[2] >= 0:
        omega_case = np.arccos((np.dot(N_vec, e_vec)/(N * e)))
    else:
        omega_case = 2 * np.pi - np.arccos((np.dot(N_vec, e_vec)/(N * e))) 
    omega_case = omega_case * 180 / np.pi #fifth orbital element
    if np.dot(N_vec, e_vec) > 0:
        print("Argument of perigee is in first or fourth quadrant")
    else:
        print("Argument of perigee is in second or third quadrant")
    #Calculating eccentric anomaly
    #EA = np.arctan2(np.sqrt((a * (1 - e**2)) / miu) * np.dot(r_vec, v_vec), (a * (1 - e**2)) - r)
    #Calculating true anomaly
    if v_rad >= 0:
        theta = np.arccos((np.dot(e_vec, r_vec)/(e * r)))
    else:
        theta = 2 * np.pi - np.arccos((np.dot(e_vec, r_vec)/(e * r))) 
    theta = theta * 180 / np.pi #sixth orbital element
    if np.dot(e_vec, r_vec) > 0:
        print("True anomaly is in first or fourth quadrant")
    else:
        print("True anomaly of perigee is in second or third quadrant")
    #Calculating mean anomaly
    #MA = EA - e * np.sin(EA)
    #MA = MA  *180 / np.pi
    #Calculationg radius of perigee and radius of apogee
    rp = (h**2 / miu) * (1 / (1 + e * np.cos(0)))
    ra = (h**2 / miu) * (1 / (1 + e * np.cos(np.pi)))
    #Calculating semi major axis
    a = 0.5 * (rp + ra)
    #Calculating mean motion
    n = np.sqrt(miu / (a**3))
    #Calculating period of an orbit
    T = (2 * np.pi) / n
    T = T / 3600
        
    return a, h, i, omega_capt, e, omega_case, theta, rp, ra, T

#testing using the r and v value from Orbital Mechanics page 161 Example 4.3
r_test = np.array([-6045, -3490, 2500])
v_test = np.array([-3.457, 6.618, 2.533])

a, h, i, omega_capt, e, omega_case, theta, rp, ra, T = sv2coe(r_test, v_test)

print("")
print("The semi major axis is ", np.absolute(a), u.km)
print("The inclination is ", i, u.degree)
print("The eccentricity is ", e)
print("The right ascension of ascending node is ", omega_capt, u.degree)
print("The argument of perigee is ", omega_case, u.degree)
print("The true anomaly is ", theta, u.degree)
print("")

print(a, h, i, omega_capt, e, omega_case, theta, rp, ra, T)

# In[ ]:




