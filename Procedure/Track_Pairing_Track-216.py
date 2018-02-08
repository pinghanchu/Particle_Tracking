#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp
import math
from scipy.misc import imsave
from scipy.misc import imread
import re
################################################
# Read frames which have been data-cleaned (wire_particle_tracking_datacleaning). 
################################################
shot = 216
picL = 18333
picR = 18158
BR = [[5.7149E-08,1.1882E-05,-2.4321E-03],[8.8831E-06,-5.3566E-07,-3.8744E-02],[-1.6469E-03,3.3774E-02,9.9867E-01]]
TR = [[-5.7704E-08,1.1998E-05,-2.1514E-03],[8.9694E-06,5.4087E-07,-3.9328E-02],[-1.7813E-03,3.3894E-02,9.9865E-01]]
TL = [[1.1040E-07,2.2954E-05,-4.1584E-03],[1.7160E-05,-1.0348E-06,6.8654E-02],[-3.4505E-03,-7.3661E-02,9.9490E-01]]
BL = [[-1.0175E-07,2.1155E-05,-4.2910E-03],[1.5815E-05,9.5369E-07,6.2907E-02],[-2.8931E-03,-6.8254E-02,9.9567E-01]]
F_10point = [[  8.34360564e-08,   1.89022714e-05,  -3.73041652e-03],[  1.94501963e-05,   1.19512678e-06,  -5.90813947e-02],[ -3.16737708e-03,   5.30589622e-02,   9.96830069e-01]]
F_big = [[  1.75345969e-07,   1.87353876e-05,  -3.69560210e-03], [  1.92508173e-05,   9.79548859e-07,  -5.72239187e-02],[ -3.15390413e-03,   5.12220927e-02,   9.97034661e-01]]
F = F_big
Matr = F
MatrixL = np.array(Matr)
MatrixR = MatrixL.transpose()
Linex = np.arange(0,384,1)
v0L = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL_sum.tif'.format(shot,shot,picL))
v0R = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR_sum.tif'.format(shot,shot,picR))
bk0L = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameL0.tif'.format(shot,shot,picL))
bk0R = imread('./Data/Shot{}/Clean_Data_Shot{}_Cam_{}/FrameR0.tif'.format(shot,shot,picR))

tL = pd.read_csv('./Data/Shot{}/trackL3_frame_inv.csv'.format(shot))
tR = pd.read_csv('./Data/Shot{}/trackR3_frame_inv.csv'.format(shot))
listL = tL['particle'].unique()
listL = listL.astype(int)
listR = tR['particle'].unique()
listR = listR.astype(int)
##Try to pair particle track between two cameras.
##Loop particle track in the left camera
for ilistL in listL:
    Pair = []
    ipL0 = int(ilistL)
    TrackL = tL[tL['particle']==ipL0]
    lenL = len(TrackL)
    if ipL0>=0 and lenL>400:
        print(ipL0,lenL)    
        xmaxL = np.max(TrackL['x'])
        xminL = np.min(TrackL['x'])
        ymaxL = np.max(TrackL['y'])
        yminL = np.min(TrackL['y'])
        for ilistR in listR:
            MassDiffSum = 0
            FFSum = 0
            FrameCount = 0
            MassDiffAve = 0
            FFAve = 0
            FFSmallCount = 0
            ipR0 = ilistR
            TrackR = tR[tR['particle']==ipR0] 
            lenR = len(TrackR)
            xmaxR = np.max(TrackR['x'])
            xminR = np.min(TrackR['x'])
            ymaxR = np.max(TrackR['y'])
            yminR = np.min(TrackR['y'])
            for iframeL in TrackL['frame'].unique():
                for iframeR in TrackR['frame'].unique():
                    if(iframeL == iframeR):
                        FrameCount = FrameCount+1
                        TrackL1=TrackL[TrackL['frame']==iframeL]
                        TrackR1=TrackR[TrackR['frame']==iframeR]
                        imassL = TrackL1['mass'].iloc[0]
                        imassR = TrackR1['mass'].iloc[0]
                        MassDiff = (imassL-imassR)*(imassL-imassR)
                        MassDiffSum = MassDiffSum+MassDiff
                        ixL = TrackL1['x'].iloc[0]
                        iyL = TrackL1['y'].iloc[0]  
                        ixR = TrackR1['x'].iloc[0]
                        iyR = TrackR1['y'].iloc[0]                         
                        ipL = [ixL,iyL,1]                       
                        ipR = [ixR,iyR,1]
                        LineL = np.dot(ipR,MatrixL) 
                        LineR = np.dot(ipL,MatrixR) 
                        #LineLy = (-LineL[2]-LineL[0]*Linex)/LineL[1]
                        #LineRy = (-LineR[2]-LineR[0]*Linex)/LineR[1]  
                        dlr = math.fabs(np.dot(ipL,LineL))/math.sqrt(LineL[0]*LineL[0]+LineL[1]*LineL[1])
                        drl = math.fabs(np.dot(ipR,LineR))/math.sqrt(LineR[0]*LineR[0]+LineR[1]*LineR[1])
                        FF = dlr + drl
                        if(FF<5):
                            FFSmallCount = FFSmallCount+1
                        FFSum = FFSum + FF
            if(FrameCount>0):
                MassDiffAve = math.sqrt(MassDiffSum/FrameCount)
                FFAve = FFSum/FrameCount
                templist = [ipL0,ipR0,lenL,lenR,xminL,xmaxL,yminL,ymaxL,xminR,xmaxR,yminR,ymaxR,FrameCount,MassDiffAve,FFAve,FFSmallCount]
                Pair.append(templist) 
                P1 = pd.DataFrame(Pair)
                P1.columns = ['PIDL','PIDR','LenL','LenR','Xmin L','Xmax L','Ymin L','Ymax L','Xmin R','Xmax R','Ymin R','Ymax R','Frame Count','Mass Diff','FF Ave','FF Small Count']
                P1.to_csv("./Data/Shot{}/TrackPair_PID_{}.csv".format(shot,ipL0))





