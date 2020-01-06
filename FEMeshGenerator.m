


function [meshOutput] = FEMeshGenerator(stlfile)
%%
% clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.3;
faceAlpha2=1;
cMap=gjet(4); 
patchColor=cMap(1,:);
markerSize=10; 

stlfile = 'STLGeometry_14_Supine_L2L3.stl';
[FV] = stlread(stlfile);
[F,V] = stlread(stlfile);

%%
% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(fileparts(filePath)),'MRI_Image_Analysis','Functions','GIBBON-master','GIBBON-master','data','temp');
modelName=fullfile(savePath,'tetgenmodel');

%% Create Outer Surface
% Specify dimensions of cube
stlfileouter = 'STLGeometryOuter_14_Supine_L2L3.stl';
[FV2] = stlread(stlfileouter);
[F2,V2] = stlread(stlfileouter);


%%
% Plotting inner disc model
cFigure; hold on;
title('Inner Disc Surface Model','FontSize',fontSize);
gpatch(F,V,patchColor,'k',faceAlpha1);
% patchNormPlot(F,V);
camlight headlight;
axisGeom(gca,fontSize); 
drawnow;

%%
% Plotting Outer disc model
cFigure; hold on;
title('Outer Cube Surface Model','FontSize',fontSize);
gpatch(F2,V2,patchColor,'k',faceAlpha1);
% patchnormPlot(F2,V2);
camlight headlight;
axisGeom(gca,fontSize);
drawnow;


%%
% Join geometries together

[F,V,C]=joinElementSets({F,F2},{V,V2});


%%
% Visualize
cFigure;
hold on;
gpatch(F,V,C,'k',0.5);
colormap(gjet(2)); icolorbar;
axisGeom;
camlight headlight;
drawnow;


%%
% DEFINE FACE BOUNDARY MARKERS
faceBoundaryMarker=ones(size(F,1),1);

%%
% Define region points
[V_regions]=getInnerPoint(F,V);

%%
% Define hole points
V_holes=[];

%% 
% Regional mesh volume parameter
[regionA]=tetVolMeanEst(F,V); %Volume for regular tets


%% 
% CREATING THE INPUT STRUCTURE
stringOpt='-pq1.2AaY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=fliplr(F);
inputStruct.Nodes=V;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA*1;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;


%% 
% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput);


end


