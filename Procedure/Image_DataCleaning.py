#!/usr/bin/env python 

# coding: utf-8

# In[1]:


from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt

# change the following to %matplotlib notebook for interactive plotting
#get_ipython().magic(u'matplotlib inline')

# Optionally, tweak styles.
mpl.rc('figure',  figsize=(10, 6))
#mpl.rc('image', cmap='gray')
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
from scipy.misc import imsave
import pims
import trackpy as tp
from pims import Frame


# In[ ]:


import trackpy
import trackpy.diag

trackpy.__version__
#trackpy.diag.dependencies()


# In[ ]:


shot=235
pic = 18333
v = pims.TiffStack('/Users/pinghanchu/Documents/Git/Data/Shot{}/Shot{}_Cam_{}.tif'.format(shot,shot,pic))
size = len(v)
v0 = tp.bandpass(v[0],0,300,threshold=5,truncate=4)
imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL0.tif".format(shot,shot,pic),v0)
zero = tp.bandpass(v[300],0,300,threshold=5,truncate=10)-v0
zero[zero>0]=0
zero[zero<0]=0
a = zero
b = zero
comb_index = 1
run_range = int(size/comb_index)
for iv1 in range(0,run_range):
    b = zero
    for iv2 in range(0,comb_index):
        iv = iv1*comb_index +iv2
        vi = tp.bandpass(v[iv],0,300,threshold=5,truncate=10)-v0
        vi[vi < 0] = 0
        vi = tp.bandpass(vi,0,300,threshold=5,truncate=10)
        vi = tp.bandpass(vi,0,300,threshold=10,truncate=10)
        b = b+vi
        a = a+vi
    #b[b>25] = 255
    imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/frame/frame_{}.tif".format(shot,shot,pic,iv1),b)
a = tp.bandpass(a,0,300,threshold=5,truncate=6)
imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL_sum.tif".format(shot,shot,pic),a)


# In[ ]:


#c = a/200+v0
#plt.imshow(c)
#imsave("/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_v0.tif".format(pic),c)
#d = a
#d[d>0]=10
#imsave("/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_{}_white.tif".format(pic,pic),d)


# In[ ]:


pic = 18158
v = pims.TiffStack('/Users/pinghanchu/Documents/Git/Data/Shot{}/Shot{}_Cam_{}.tif'.format(shot,shot,pic))
v0 = tp.bandpass(v[0],0,300,threshold=5,truncate=4)
imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR0.tif".format(shot,shot,pic),v0)
#plt.imshow(v0)
zero = tp.bandpass(v[300],0,300,threshold=5,truncate=10)-v0
zero[zero>0]=0
zero[zero<0]=0
a = zero
b = zero
for iv1 in range(0,run_range):
    b = zero
    for iv2 in range(0,comb_index):
        iv = iv1*comb_index+iv2
        vi = tp.bandpass(v[iv],0,300,threshold=5,truncate=10)-v0
        vi[vi < 0] = 0
        vi = tp.bandpass(vi,0,300,threshold=5,truncate=10)
        vi = tp.bandpass(vi,0,300,threshold=10,truncate=10)
        b = b+vi
        a = a+vi
    #b[b>25] = 255
    imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/frame/frame_{}.tif".format(shot,shot,pic,iv1),b)
a = tp.bandpass(a,0,300,threshold=5,truncate=6)
imsave("/Users/pinghanchu/Documents/Git/Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR_sum.tif".format(shot,shot,pic),a)
#c = a/200+v0
#plt.imshow(c)
#imsave("/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_v0.tif".format(pic),c)
#d = a
#d[d>0]=10
#imsave("/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_{}_white.tif".format(pic,pic),d)

