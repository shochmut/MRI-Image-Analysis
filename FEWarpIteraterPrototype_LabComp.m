%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               %%%%%%%%%%%%   %%%%%%%%%%%%%    %%%%%%%%%%%%%%            %
%               %              %                %           %             %
%               %              %                %          %              %
%               %%%%%%%%%      %%%%%%%%%%       %%%%%%%%%%%     io        %
%               %              %                %          %              %
%               %              %                %           %             %
%               %              %%%%%%%%%%%%%    %%%%%%%%%%%%%%            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skyler Hochmuth
% Colorado State University 
% Walter Scott Junior College of Engineering
function [] = FEWarpIteraterPrototype_LabComp(E_youngs,v_poisson)

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting force
febioLogFileName_warptemplate=[febioFebFileNamePart,'_warptemplate_out.txt']; %Log file name for exporting stress
febioLogFileName_warptarget=[febioFebFileNamePart,'_warptarget_out.txt']; %Log file name for exporting stiffness
febioLogFileName_warpenergy=[febioFebFileNamePart,'_warpenergy_out.txt'];
febioLogFileName_warpforce=[febioFebFileNamePart,'_warpforce_out.txt'];


%Specifying dimensions and number of elements
cubeSize=50; 
sampleWidth=50; %Width 
sampleThickness=50; %Thickness 
sampleHeight=20; %Height
pointSpacings=2*ones(1,3); %Desired point spacing between nodes
numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3

%Define the material density, E (young's modulus), and v (poisson ratio)
density = 1;
E_youngs = E_youngs;
v_poisson = v_poisson;

%Define applied displacement 
appliedStrain=0.4; %Linear strain (Only used to compute applied stretch)
loadingOption='tension'; % or 'tension'
switch loadingOption
    case 'compression'
        stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
    case 'tension'
        stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
end
displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight; %The displacement magnitude

%Material parameter set
c1=0.7; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=500; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
formulationType='uncoupled'; %coupled

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%% Run tetgen meshing function to get inner meshed object
[meshOutput] = FEMeshGenerator();

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

meshStruct = meshOutput;



%Access elements, nodes, and faces from the structure
E=meshStruct.elements; %The elements 
V=meshStruct.nodes; %The nodes (vertices)
Fb=meshStruct.facesBoundary; %The boundary faces
Cb=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
elementMaterialIndices=ones(size(E,1),1); %Element material indices

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; 

subplot(1,2,1); hold on; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',faceAlpha1); 
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);

drawnow;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=fewarpStructTemplate_LabComp;

%febio_spec version 
febio_spec.ATTR.version='2.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis.ATTR.type='static';
febio_spec.Control.title='Cube analysis';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.dtol=0.001;
febio_spec.Control.etol=0.01;
febio_spec.Control.rtol=1;
febio_spec.Control.lstol=0.9;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;



febio_spec.Globals.Constants.T=0;
febio_spec.Globals.Constants.R=0;
febio_spec.Globals.Constants.Fc=0;

%Material section
switch formulationType
    case 'coupled'
        febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.density=density;
        febio_spec.Material.material{1}.E=E_youngs;
        febio_spec.Material.material{1}.v=v_poisson;
    case 'uncoupled'
        febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.density=density;
        febio_spec.Material.material{1}.E=E_youngs;
        febio_spec.Material.material{1}.v=v_poisson;
end


%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Nucleus_Pulposus'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;


%% Output section 
% -> log file
% febio_spec.Output.logfile.ATTR.file=febioLogFileName;
% febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
% febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
% febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
% febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);
% 
% febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
% febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
% febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);
% 
% febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_warptemplate;
% febio_spec.Output.logfile.element_data{2}.ATTR.data='cxxxx;cxxyy;cyyyy;cxxzz;cyyzz;czzzz;cxxxy;cyyxy;czzxy;cxyxy;cxxyz;cyyyz;czzyz;cxyyz;cyzyz;cxxxz;cyyxz;czzxz;cxyxz;cyzxz;cxzxz';
% febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{2}.VAL=1:size(E,1);
% 
% febio_spec.Output.logfile.element_data{3}.ATTR.file=febioLogFileName_warptarget;
% febio_spec.Output.logfile.element_data{3}.ATTR.data='cxxxx;cxxyy;cyyyy;cxxzz;cyyzz;czzzz;cxxxy;cyyxy;czzxy;cxyxy;cxxyz;cyyyz;czzyz;cxyyz;cyzyz;cxxxz;cyyxz;czzxz;cxyxz;cyzxz;cxzxz';
% febio_spec.Output.logfile.element_data{3}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{3}.VAL=1:size(E,1);
% 
% febio_spec.Output.logfile.element_data{4}.ATTR.file=febioLogFileName_warpenergy;
% febio_spec.Output.logfile.element_data{4}.ATTR.data='cxxxx;cxxyy;cyyyy;cxxzz;cyyzz;czzzz;cxxxy;cyyxy;czzxy;cxyxy;cxxyz;cyyyz;czzyz;cxyyz;cyzyz;cxxxz;cyyxz;czzxz;cxyxz;cyzxz;cxzxz';
% febio_spec.Output.logfile.element_data{4}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{4}.VAL=1:size(E,1);
% 
% febio_spec.Output.logfile.element_data{5}.ATTR.file=febioLogFileName_warpforce;
% febio_spec.Output.logfile.element_data{5}.ATTR.data='cxxxx;cxxyy;cyyyy;cxxzz;cyyzz;czzzz;cxxxy;cyyxy;czzxy;cxyxy;cxxyz;cyyyz;czzyz;cxyyz;cyzyz;cxxxz;cyyxz;czzxz;cxyxz;cyzxz;cxzxz';
% febio_spec.Output.logfile.element_data{5}.ATTR.delim=',';
% febio_spec.Output.logfile.element_data{5}.VAL=1:size(E,1);



%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='internal';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [~, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp));
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    % Importing element stress from a log file
    [time_mat, E_stress_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress)); 
    time_mat=[0; time_mat(:)]; %Time
    stress_cauchy_sim=[0; mean(squeeze(E_stress_mat(:,end,:)),1)'];
    
    %% Importing stiffness data
    % Importing element stiffness tensors from a log file
%     [~,stiffness_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_stiffness)); 
%     
%     stiffness_mat=stiffness_mat(:,2:end,end); %Final stiffness state
%     
%     stiffness_mat_voigt=stiffness_mat(:,[1  2  4  11 16 7;...
%                                          2  3  5  12 17 8;...
%                                          4  5  6  13 18 9;...
%                                          11 12 13 15 20 14;...
%                                          16 17 18 20 21 19;...
%                                          7  8  9  14 19 10]);
%     stiffness_mat_voigt=reshape(stiffness_mat_voigt',6,6,size(stiffness_mat_voigt,1));
%     stiffness_mat_voigt=reshape(mat2cell(stiffness_mat_voigt,6,6,...
%         ones(size(stiffness_mat_voigt,3),1)),[size(stiffness_mat,1),1]);
%     
%     stiffness_mat_kelvin=stiffness_mat_voigt; 
%     for q=1:1:numel(stiffness_mat_voigt)
%         cVoigt=stiffness_mat_voigt{q};
%         c=voigtUnMap(cVoigt);
%         cKelvin=kelvinMap(c);
%         stiffness_mat_kelvin{q}=cKelvin;
%     end
     
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp=gpatch(Fb,V_def,CF,'k',1); %Add graphics object to animate
    gpatch(Fb,V,0.5*ones(1,3),'k',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis([min(V_def(:,1)) max(V_def(:,1)) min(V_def(:,2)) max(V_def(:,2)) min(V_def(:,3)) max(V_def(:,3))]); %Set axis limits statically
    view(130,25); %Set view direction
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        [CF]=vertexToFaceMeasure(Fb,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    
    %% 
    % Calculate the simulated applied uniaxial stretch
    
    DZ_set=N_disp_mat(bcPrescribeList,end,:); %Z displacements of the prescribed set
    DZ_set=mean(DZ_set,1); %Calculate mean Z displacements across nodes
    stretch_sim=(DZ_set(:)+sampleHeight)./sampleHeight; %Derive stretch
        
    %%    
    % Visualize stress-stretch curve
    
    cFigure;
    hold on;    
    title('Uniaxial stress-stretch curve','FontSize',fontSize);
    xlabel('\lambda Stretch [.]','FontSize',fontSize); ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize); 
    
    plot(stretch_sim(:),stress_cauchy_sim(:),'r-','lineWidth',lineWidth);
    
    view(2); axis tight;  grid on; axis square; box on; 
    set(gca,'FontSize',fontSize);
    drawnow;
    
end
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
