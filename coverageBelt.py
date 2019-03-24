#!/usr/bin/env python
# coding: utf-8

# In[4]:


import twoBodyTool
from scipy.integrate import odeint
import numpy as np
from twoBodyTool import F, F1, sv2coe, RE, J2, miu, rate, coverageBelt

A = np.array([600, 700, 800, 900, 1000])
B = np.array([8, 8, 8, 8, 8])

slantRange, nadirAngle, centralAngle, coverageArea, coverPercent = coverageBelt(A, B)
