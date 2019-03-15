#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np

#G_const for gravitational constant, ME for mass of earth, RE for radius of earth
G_const = 6.67408 * 10**(-11)
ME = 5.972364730419773 * 10**24
RE = 6378100
#Declaring value of miu
miu = (G_const * ME) * 10**(-9)
#Declaring teh value of perturbation J2
J2 = 1.082629 * 10**(-3)

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
#   a = the magnitude of semimajor axis (in km)
#   i = the magnitude of inclination (in degree)
#   omega_capt = the magnitude of right ascencion of the ascending node (in degree)
#   e = the magnitude of eccentricity (no unit)
#   omega_case = the magnitude of argument of perigee (in degree)
#   theta = the magnitude of true anomaly (in degree)
#   h = the magnitude of specific angular momentum (in km^2/s)
#   T = the magnitude of period the orbit (in hour)
# all the output will be in form of output array
# The output of the function needs to be recalled as an array function of 
# [lists of desired output] = sv2coe(position vetor, velocity vector)
def sv2coe(r_vec, v_vec):
    #Calculating distance from r_vec
    r = np.sqrt(np.dot(r_vec, r_vec))
    #Calculating speed from v_vec
    v = np.sqrt(np.dot(v_vec, v_vec))
    #Calculating radial velocity
    v_rad = (np.dot(r_vec, v_vec)) / r
    #if v_rad > 0:
    #    print("Satellite is flying away from perigee")
    #else:
    #    print("Satellite is flying towards perigee")
    #Calculating specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    #Calculating magnitude of specific angular momentum
    h = np.sqrt(np.dot(h_vec, h_vec)) #first orbital element
    #Calculating inclination
    i = np.arccos(h_vec[2] / h) 
    i = i*180/np.pi #second orbital element
    #if i <= 90:
    #    print("The orbit is retrograde")
    #elif 90 < i < 180:
    #    print("The orbit is posigrade")
    #else:
    #    print("The orbit is greater than 180 and in quadrant ambiguity")
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
    #if np.dot(N_vec, e_vec) > 0:
    #    print("Argument of perigee is in first or fourth quadrant")
    #else:
    #    print("Argument of perigee is in second or third quadrant")
    #Calculating eccentric anomaly
    #EA = np.arctan2(np.sqrt((a * (1 - e**2)) / miu) * np.dot(r_vec, v_vec), (a * (1 - e**2)) - r)
    #Calculating true anomaly
    if v_rad >= 0:
        theta = np.arccos((np.dot(e_vec, r_vec)/(e * r)))
    else:
        theta = 2 * np.pi - np.arccos((np.dot(e_vec, r_vec)/(e * r))) 
    theta = theta * 180 / np.pi #sixth orbital element
    #if np.dot(e_vec, r_vec) > 0:
    #    print("True anomaly is in first or fourth quadrant")
    #else:
    #    print("True anomaly of perigee is in second or third quadrant")
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
    #Calculating semi-latus rectum
    p = a * (1 - e**2)    
    return [a, i, omega_capt, e, omega_case, theta, h, T, n, p]

# defining function coe2sv(h, i, omega_capt, e, omega_case, theta) for transforming state vector to classical orbital elements
# in the following definition of function sv2coe, sv is a state vector and coe is classical orbital elements
# the input element of the function sv2coe that need to be defined consists of:
#   h = the magnitude of specific angular momentum (in km^2/s)
#   i = the magnitude of inclination (in degree)
#   omega_capt = the magnitude of right ascencion of the ascending node (in degree)
#   e = the magnitude of eccentricity (no unit)
#   omega_case = the magnitude of argument of perigee (in degree)
#   theta = the magnitude of true anomaly (in degree)
#   all the inputs should be in each correct unit
# 
# the possible output of this function are as follows:
#   r_peri = the position vector in perifocal coordinates (in km)
#   v_peri = the velocity vector in perifocal coordinates (in km/s)
#   Q_geo = the matrix Q_geo transformation from perifocal to geocentric equatorial coordinates (no unit)
#   r_geo = the position vector in geocentric frame (in km)
#   v_geo = the velocity vector in geocentric frame (in km/s)
# all the output will be in form of output array
# The output of the function needs to be recalled as an array function of 
# [lists of desired output] = coe2sv(h, i, omega_capt, e, omega_case, theta)
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
    
    return [r_peri, v_peri, Q_geo, r_geo, v_geo]

from scipy.integrate import odeint

# In the following definition of orbit propagator derivative function F,
# s is a state vector whose components represent the following:
#   s[0] = Horizontal or x position
#   s[1] = Vertical or y position
#   s[2] = Lateral or z position
#   s[3] = Horizontal velocity
#   s[4] = Vertical velocity
#   s[5] = Lateral velocity
# while t is the period of the orbit (in second) that needs to be defined more further; t = np.linspace(start,end,interval)
#
# The input for the defined derivative function consists of:
#   x0 = Initial position in x axis
#   y0 = Initial position in y axis
#   z0 = Initial position in z axis
#   vx0 = Initial velocity in x axis
#   vy0 = Initial velocity in y axis
#   vz0 = Initial velocity in z axis
# The solution of the derivative function will be using odeint() function, where will be declared as:
# solution = odeint(F,[x0,y0,z0,vx0,vy0,vz0],t)
def F(s,t):
    a = -miu*s[0]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2)
    b = -miu*s[1]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2)
    c = -miu*s[2]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2)
    return [s[3],s[4], s[5], a, b, c]

# In the following definition of orbit propagator derivative function F1 or Perturbation Accelerations due to J2,
# s is a state vector whose components represent the following:
#   s[0] = Horizontal or x position
#   s[1] = Vertical or y position
#   s[2] = Lateral or z position
#   s[3] = Horizontal velocity
#   s[4] = Vertical velocity
#   s[5] = Lateral velocity
# while t is the period of the orbit (in second) that needs to be defined more further; t = np.linspace(start,end,interval)
#
# The input for the defined derivative function consists of:
#   x0 = Initial position in x axis
#   y0 = Initial position in y axis
#   z0 = Initial position in z axis
#   vx0 = Initial velocity in x axis
#   vy0 = Initial velocity in y axis
#   vz0 = Initial velocity in z axis
# The solution of the derivative function will be using odeint() function, where will be declared as:
# solution = odeint(F,[x0,y0,z0,vx0,vy0,vz0],t)
#
# Plotting the orbit from orbit propagator
# 1. Define the period (t) of the orbit in second
# 2. Declare the solution function using odeint() function
# 3. Declare a void array, e.g. results = []
# 4. Example of getting the array of integral results:
#    for data in range(0,1000):
#      r_test = np.squeeze(solution[data:data+1,0:3])
#      v_test = np.squeeze(solution[data:data+1,3:6])
#      t1 = np.linspace(0,8200,data+1)
#      results.append(sv2coe(r_test, v_test))
#    results = np.array(results)
# 5. Example of command plotting:
#    import matplotlib.pyplot as plt
#    x1 = results[:,0]
#    y1 = t1
#    plt.plot(y1,x1, '.')
#    plt.axis('equal')
#    plt.show()
def F1(s,t):
    Pt = (3/2)*(miu*J2*((RE*10**(-3))**2))
    r = (s[0]**2 + s[1]**2 + s[2]**2)**(1/2)
    a = -miu*s[0]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2) - Pt * ((1 - 5*(s[2]**2/r**2))*(s[0]/r**5))
    b = -miu*s[1]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2) - Pt * ((1 - 5*(s[2]**2/r**2))*(s[1]/r**5))
    c = -miu*s[2]/(s[0]**2 + s[1]**2 + s[2]**2)**(3/2) - Pt * ((3 - 5*(s[2]**2/r**2))*(s[2]/r**5))
    return [s[3],s[4], s[5], a, b, c]

# The rate() function is defining the value of the mean classical orbital elements rate of change,
# the input components of the function is described as following:
# i = inclination of the orbit (in degree)
# e = eccentricity of the orbit (no unit)
# n = mean motion of the orbit (no unit)
# p = semi-latus rectum/orbit parameter (no unit)
#
# The output components of the function comprises of:
# semiMajRate = rate of semi major axis (in km/second)
# eccenRate = rate of eccentricity (no unit)
# incliRate = rate of inclination (in degree/second)
# omgRate = rate of right ascencion of the ascending node (in degree/second)
# omgCaseRate = rate of argument of perigee (in degree/second)
# meanAnomaly = rate of mean anomaly (in degree/second)
def rate(i, e, n, p):
    rateFactor = n*J2*((RE*10**(-3)/p)**2)
    semiMajRate = 0
    eccenRate = 0
    incliRate = 0
    omgRate = -(3/2)*(rateFactor*np.cos(i*np.pi/180))
    omgCaseRate = (3/4)*rateFactor*(5*(np.cos(i*np.pi/180)**2)-1)
    meanAnomaly = (n + (3/4)*np.sqrt(1-e**2)*rateFactor*(3*(np.cos(i*np.pi/180)**2)-1))
    return semiMajRate, eccenRate, incliRate, omgRate, omgCaseRate, meanAnomaly

# In[ ]:




