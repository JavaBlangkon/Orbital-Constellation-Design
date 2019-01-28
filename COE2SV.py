#!/usr/bin/env python
# coding: utf-8

# In[23]:


#astropy is a python library for astronomy
#G for gravitational constant, M_earth for mass of earth, R_earth for radius of earth
from astropy.constants import G, M_earth, R_earth
from astropy import units as u
import numpy as np

#declaring value of miu and radius of earth
miu = (G.value*M_earth.value)*10**(-9)
Re = R_earth.value

#declaring unit vector
i_unit = np.array([1, 0, 0])
j_unit = np.array([0, 1, 0])
k_unit = np.array([0, 0, 1])

#defining function coe2sv for transforming classical orbital elements to state vector
def coe2sv(r_vec,v_vec):
    #Calculating distance from r_vec
    r = np.sqrt(np.dot(r_vec, r_vec))
    #Calculating speed from v_vec
    v = np.sqrt(np.dot(v_vec, v_vec))
    #Calculating radial velocity
    v_rad = (np.dot(r_vec, v_vec))/r
    if v_rad > 0:
        print("Satellite is flying away from perigee")
    else:
        print("Satellite is flying towards perigee")
    #Calculating specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    #Calculating magnitude of specific angular momentum
    h = np.sqrt(np.dot(h_vec, h_vec)) #first orbital element
    #Calculating inclination
    i = np.arccos(h_vec[2]/h) 
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
        omega_capt = np.arccos(N_vec[0]/N)
    else:
        omega_capt = 2*np.pi - np.arccos(N_vec[0]/N) 
    omega_capt = omega_capt*180/np.pi #third orbital element
    #Calculating eccentricity vector
    e_vec = (1/miu)*(np.dot((v**2 - (miu/r)), r_vec) - (np.dot((r*v_rad), v_vec)))
    #Calculating magnitude of eccentricity
    e = np.sqrt(np.dot(e_vec, e_vec)) #fourth orbital element
    #Calculating semi major axis
    a = -miu/(2*e)
    #Calculating argument of perigee
    if e_vec[2] >= 0:
        omega_case = np.arccos((np.dot(N_vec, e_vec)/(N*e)))
    else:
        omega_case = 2*np.pi - np.arccos((np.dot(N_vec, e_vec)/(N*e))) 
    omega_case = omega_case*180/np.pi #fifth orbital element
    if np.dot(N_vec, e_vec) > 0:
        print("Argument of perigee is in first or fourth quadrant")
    else:
        print("Argument of perigee is in second or third quadrant")
    #Calculating true anomaly
    if v_rad >= 0:
        theta = np.arccos((np.dot(e_vec, r_vec)/(e*r)))
    else:
        theta = 2*np.pi - np.arccos((np.dot(e_vec, r_vec)/(e*r))) 
    theta = theta*180/np.pi #sixth orbital element
    if np.dot(e_vec, r_vec) > 0:
        print("True anomaly is in first or fourth quadrant")
    else:
        print("True anomaly of perigee is in second or third quadrant")
        
    return a, h, i, omega_capt, e, omega_case, theta

#testing using the r and v value from Orbital Mechanics page 161 Example 4.3
r_test = np.array([-6045, -3490, 2500])
v_test = np.array([-3.457, 6.618, 2.533])

a, h, i, omega_capt, e, omega_case, theta = coe2sv(r_test, v_test)

print("")
print("The semi major axis is ", np.absolute(a), u.km)
print("The inclination is ", i, u.degree)
print("The eccentricity is ", e)
print("The right ascension of ascending node is ", omega_capt, u.degree)
print("The argument of perigee is ", omega_case, u.degree)
print("The true anomaly is ", theta, u.degree)
print("")

print(a, h, i, omega_capt, e, omega_case, theta)
print(miu) #need more explanation on this, since the value in the book is different from real calculation assumption


# In[ ]:




