

clear
clc
close all

%% First add the subpath of the directory: (You must be at the root of the
% directory)
addpath(genpath('.'));


ImgStack = uigetfile('../.tif','Select initial image .tiff file');


t = Tiff(ImgStack,'r');                             % open and read tif file

%% Convert the image stack into a 3d .raw file for input into FEBio
fid = fopen(strcat('raw','_',ImgStack,'.raw'),'wt');
finalraw = fwrite(fid,ImgStack,'uint8');
fclose(fid);

