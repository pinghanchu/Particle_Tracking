#!/usr/bin/env python                                                           

# coding: utf-8

# In[1]:


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp
import trackpy.predict
from scipy.misc import imsave
from scipy.misc import imread
import math
################################################
# Read frames which have been data-cleaned (wire_particle_tracking_datacleaning). 
################################################
shot = 232
picL = 18333
picR = 18158
FrameL = pims.ImageSequence('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/invframe/frame_*.tif'.format(shot,shot,picL))
FrameR = pims.ImageSequence('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/invframe/frame_*.tif'.format(shot,shot,picR))
#FrameL = pims.ImageSequence('/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_white_*.tif'.format(picL))
#FrameR = pims.ImageSequence('/Users/pinghanchu/Documents/Git/Data/Clean_Data_Shot119_Cam_{}/frame_white_*.tif'.format(picR))
v0L = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL_sum.tif'.format(shot,shot,picL))
v0R = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR_sum.tif'.format(shot,shot,picR))
bk0L = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL0.tif'.format(shot,shot,picL))
bk0R = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR0.tif'.format(shot,shot,picR))


# In[ ]:


###################################
#Locate Features; single frame test
###################################
#init_index = 100
#f = tp.locate(FrameL[init_index], 3, minmass=5) 
#Show points located. It will be better to have more points rather than missing points.
#plt.figure(figsize=[12,12])  # make a new figure
#tp.annotate(f, FrameL[init_index]);


# In[ ]:


#Show 'mass' distribution
#fig, ax = plt.subplots()
#ax.hist(f['mass'], bins=20)
# Optionally, label the axes.
#ax.set(xlabel='mass', ylabel='count');
#plt.show()


# In[ ]:


######################################################################
#Locate Features; apply the same parameter to all frames
######################################################################
# Left frame
pred = trackpy.predict.NearestVelocityPredict()
#pred = trackpy.predict.ChannelPredict(0.5, 'x', minsamples=3)
fL = tp.batch(FrameL, 3, minmass=5);
tL = pred.link_df(fL, 3, memory=11,  diagnostics=True)
tL.to_csv('./Data/Shot{}/trackL_frame_inv.csv'.format(shot))


# In[ ]:


#tL = pd.read_csv('./Data/Shot{}/trackL_frame_inv.csv'.formate(shot))
tL.head()


# In[ ]:


plt.figure(figsize=[12,12])
tp.plot_traj(tL);


# In[ ]:


plt.figure(figsize=[12,12])
plt.imshow(v0L+bk0L)
#plt.scatter(tL['x'],tL['y'],s=0.3,c='g')
plt.scatter(tL['x'],tL['y'],s=0.3,c=tL['mass'])
plt.show()


# In[ ]:


# Remove tracks too few points (less than 500)
tL1 = tp.filter_stubs(tL,200)
plt.figure(figsize=[12,12])
tp.plot_traj(tL1);


# In[ ]:


plt.figure(figsize=[12,12])
plt.imshow(v0L)
plt.scatter(tL1['x'],tL1['y'],s=0.3,c=tL1['mass'])
#plt.scatter(tL1['x'],tL1['y'],s=0.3,c='r')
plt.show()


# In[ ]:


tL2 = tL1
range_limit = 10
Range = {}
for ii in tL2['particle'].unique():
    Track = tL2[tL2['particle']==ii]
    xmax = np.max(Track['x'])
    xmin = np.min(Track['x'])
    ymax = np.max(Track['y'])
    ymin = np.min(Track['y'])
    ra = math.sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin))
    Range[ii] = ra


# In[ ]:


RANGE = []
for ii in range(len(tL2)):
    RANGE.append(Range[tL2['particle'].iloc[ii]])
    #print(ii,Range[tL2['particle'].iloc[ii]])


# In[ ]:


tL2['range'] = RANGE


# In[ ]:


tL2.head()


# In[ ]:


tL3 = tL2[tL2['range']>20]


# In[ ]:


plt.figure(figsize=[12,12])
plt.imshow(v0L+bk0L)
plt.scatter(tL3['x'],tL3['y'],s=0.3,c=tL3['mass'])
#plt.scatter(tL3['x'],tL3['y'],s=0.1,c='g')
plt.savefig('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/sumTrackL.tif'.format(shot,shot,picL))
plt.show()


# In[ ]:


'''
for ii in tL3['particle'].unique():
#for ii in range(3,4):  
    #print(ii)
    Track = tL3[tL3['particle']==ii]
    xmax = np.max(Track['x'])
    xmin = np.min(Track['x'])
    ymax = np.max(Track['y'])
    ymin = np.min(Track['y'])
    plt.figure(figsize=[12,12])
    plt.imshow(v0L)
    plt.scatter(Track['x'],Track['y'],s=0.7,c=Track['mass'])
    plt.ylim(ymin-10,ymax+10)
    plt.xlim(xmin-10,xmax+10)
    plt.savefig("./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/trackL_framebk_{}.tif".format(shot,shot,picL,int(ii)))
'''


# In[ ]:


plt.figure(figsize=[12,12])
tp.plot_traj(tL3);


# In[ ]:


tL3.to_csv('./Data/Shot{}/trackL3_frame_inv.csv'.format(shot))


# In[ ]:


tL3.head()

