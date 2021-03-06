#Particle_Tracking

Programs are used to identify particle tracks.

Matlab Particle Tracking
http://site.physics.georgetown.edu/matlab/

TrackPy
https://soft-matter.github.io/trackpy/v0.3.2/

OpenPTV
http://openptv-python.readthedocs.io/en/latest/index.html

OpenCV
https://opencv.org/

Procedure:

1. Data cleaning:
./Procedure/Image_DataCleaning

It subtracts the background from each image. The background is the first image while no signal yet.
It resets the pixel to zero if the value is negative or smaller than 10. All images are stored in
../../Data/Clean_Data_Shot119_Cam_18*/frame/

It also generates the background images:
./Procedure/FrameL0.tif
./Procedure/FrameR0.tif

It also generates the combination of all images:
./Procedure/FrameL_sum.tif
./Procedure/FrameR_sum.tif

It also generates the combination of all images and make all pixel values uniform:
./Procedure/FrameL_whilte.tif
./Procedure/FrameR_whilte.tif

2. Reverse the frame sequence:
./Procedure/copy.py
Use copy.py to reverse the sequency of the images.All images are stored in 
../../Data/Clean_Data_Shot119_Cam_18*/invframe/

3. Calibration fundamental matrix:
./Procedure/Calibration_10points.ipynb
Generate calibration fundamental matrix from a few known points. 
This gives us a basic matrix for particle pairing later. 
The key is to have points spread all area of the frame as much as possible. 

4. Particle Pairing:
./Procedure/Particle_Pairing.ipynb
Use TrackPy to identify points in a frame.
Apply the fundamenatal matrix in the step 3 to identify particle pairs in the right and left frames.
Apply this on 10-ish frames to collect ~1000 points (at least).

5. Calibration fundamental matrix using big data:
./Procedure/Calibration_Npoints.ipynb
Generate calibration fundamental matrix from the points collected in the step 4.
Have the points spread all area of the frame.
This will minimize the error of the calibration fundamental matrix.

6. Tracking:
./Procedure/TrackingL.ipynb
./Procedure/TrackingR.ipynb
Use TrackPy to identify "track" from "../../Data/Clean_Data_Shot119_Cam_18*/invframe/"
using "NearestVelocityPredict()". 
Remove:	i) the tracks too few points (> 500 points); 
	ii) the tracks too short (maximum (x,y) to minimum (x,y) > 20);

Save the information to
./Procedure/trackL3_frame_inv.csv
./Procedure/trackR3_frame_inv.csv

7. Track Pairing:
./Procedure/Track_Pairing_Track.ipynb

a. Start with a track in the left frame. (PIDL) Can apply additional cut to remove short tracks.
b. Compare this track with all other tracks in the right frame
Calculate:
i. Average FF: for overlap region, calculate the average FF
ii. Average Mass Difference: for overlap region, calculate the average mass difference = sqrt(Sum^N (massL-massR)^2/N)
iii. Overlap frame number
iv. Frame overlap number with FF<5 
v. Track size(xmin,xmax,ymin,ymax)
vi. Track frame number

It generates TrackPair_PID_[i].csv for track [i] in the left camera.
It contains the corresponding parameters for PIDL and all other PIDR.

./Procedure/MergeCsvFile.ipynb
It combines all TrackPair_PID_*.csv to TrackPair_PID.csv

./Procedure/Track_Pairing_Cut.ipynb
Apply a cut on TrackPair_PID.csv to select pairs. Require i)<5, ii)<50, iii)>1000, iv>100.
It generates PairList.csv which contains the possible pair (PIDL, PIDR).
It also generates Pair_[i].csv for the track [i] in the left camera and all other pairs [j] in the right camera.

./Procedure/TrackPair_PID_Plot.ipynb
Draw four figures for one pair: Left track, right track, mass and the FF from Pair_[i].csv
This provides additional human-interactive cut.  

8. Plot and fit 3D
./Procedure/3DTrack_Fitting.ipynb
Fit 3D tracks. 