% Example of how to extract tensor statistcis form an ROI defined from a set of coordinates as those of the cluster.
%
% Script wrote fro Tanya Glozman.
%
% Copyright 2015 Franco Pestilli Indiana University

% Define an ROI from the coordinates of the cluster you have precomputed. These coordinates are expected to be in ACPC.
clusterCoors = [x,y,z];% Coordinates in ACPC
roi = dtiNewRoi('clusterNum1', 'r', clusterCoords)
 
% Load the dt6 file computed with dtiInit.m
dt = dtiLoadDt6('full/path/to/dt6.mat');

% Extract the FA, MD etc values from the diffusion image. For this you will need to have the data preprocesed with dtiInit.
[myVals1,myVals2, myVals3, myVals4, myVals5, myVals6, myVals7] = dtiGetValFromTensors(dt.dt6, roi.coords, inv(dt.xformToAcpc), 'FA');
