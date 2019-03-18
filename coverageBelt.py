#!/usr/bin/env python
# coding: utf-8

# In[4]:


import twoBodyTool
from scipy.integrate import odeint
from twoBodyTool import F, F1, sv2coe, RE, J2, miu, rate, coverageBelt

result = coverageBelt(600, 2)


# In[8]:


result*2

