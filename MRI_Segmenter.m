% Skyler Hochmuth
% Colorado State University
% Walter Scott School of Biomedical Engineering
% Spring 2019

function [final_disc, final_overlay, stacklength] = MRI_Segmenter(ImgStack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab function converts the pixelated MRI images into a coordinate%
% system. Images are then analyzed using a seed growing algorithm to find % 
% the pixels that contain the nucleus pulposus (white) contained within   %
% the annulus fibrosis (black)                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First add the subpath of the directory: (You must be at the root of the
% directory)
% addpath(genpath('.'));
t = Tiff(ImgStack,'r');

%% Run a for loop to determine the number of image slices in the tif file
for i = 1:50                                        % set up for loop 
    setDirectory(t,i)                               % set up 1st IFD
    if lastDirectory(t) == 0                        % if IFD isn't last IFD in stack force next iteration
        continue
    end
    if lastDirectory(t) == 1                        % check if it is the last IFD in stack
        stacklength = i;                           % set length of stack
        break
    end
end


%% Set up matrix of images
for i = 1:stacklength                               % initiate for loop
    A(:,:,i) = imread(ImgStack,i);
end


%% Display pre-processing of images
for i = 1:stacklength
    A_sharp1(:,:,i) = imsharpen(A(:,:,i),'Amount',3,'Radius',3);
    A_bilatfilt(:,:,i) = imbilatfilt(A_sharp1(:,:,i));
    A_histeq(:,:,i) = histeq(A_sharp1(:,:,i));
    sharp2(:,:,i) = imsharpen(A_histeq(:,:,i),'Amount',3,'Radius',8);
end
% figure
% montage({A(:,:,13),A_bilatfilt(:,:,13),A_sharp1(:,:,13),A_histeq(:,:,13),sharp2(:,:,13)},'Size',[1 5])   % display the 3 images    


%% Use a GUI to determine disc border slices
figure
x = []; y = [];                                     % preallocate empty points
montage(A,'size',[1,stacklength],'BorderSize',...
        [2,2],'BackgroundColor','r');               % set up a montage of the images
text(0,-40,'Single click the image slices where disc begins and ends');
text(0,-20,'Select the disc begin first and disc end last');
q = getframe;                                       
Q = q.cdata;
Q = frame2im(q);                                    % getframe of image 
Q = rgb2gray(Q);

[x,y] = ginput(2);  
discbegin = 1 + floor(x(1)/(length(Q(1,:))/stacklength));
discend = 1 + floor(x(2)/(length(Q(1,:))/stacklength));
close all


%% Use a GUI to determine the image that should be used for initial seeding
figure
x = []; y = [];                                     % preallocate empty points
montage(A,'size',[1,stacklength],'BorderSize',...
        [2,2],'BackgroundColor','r');               % set up a montage of the images
text(0,-10,'Double click the image slice with the most well defined disc');
q = getframe;                                       
Q = q.cdata;
Q = frame2im(q);                                    % getframe of image 
Q = rgb2gray(Q);

[x,y] = getpts;  
ImgNum = 1 + floor(x/(length(Q(1,:))/stacklength));
close all


%% GUI: Obtain the seed point for initial image
InitImg = A(:,:,ImgNum);
figure
x = []; y = [];                                     % preallocate empty points
imshow(InitImg);
text(0,-10,'Double click on seed point')
text(0,-5,'If no visible seed point press enter')
[x,y] = getpts;                             % get seed point
x = round(x);                                       % round seed point
y = round(y);
close all                                           % close images


%% Set up initial tolerance for seed algorithm
tol = 0.4;


%% Run the seed growing algorithm on initial image
[final, colored, pouthisteq] = Image_Analysis_Fxn(InitImg,x,y,tol); % run image analysis fxn and store final images
    MaskInit(:,:) = final(:,:);                       % store inital image mask
% figure
% imshow(colored)
% title('Initial Mask Generation')


%% If initial seed grow touches border then sharpen image even more
% This step only acts when the condition is met in order to prevent the
% rest of the image analysis from being innaccurate by narrowing the
% results. This step takes away some accuracy in favor of a more robust code
InitClearBord = imclearborder(MaskInit);             % if disc touches border clear
Inittrue = find(InitClearBord);                      % find images that have been cleared
    if isempty(Inittrue) == 1
        close all
        InitImg = imsharpen(InitImg,'Amount',2,'Radius',2);
        [final, colored, pouthisteq] = Image_Analysis_Fxn(InitImg,x,y,tol); % run image analysis fxn and store final images
        MaskInit(:,:) = final(:,:);                       % store inital image mask
        figure
        imshow(colored)
        title('Initial Mask Generation')
    end 

    
%% Use the obtained initial mask to find seed points in all other image slices
% This step finds the brightest pixel closest to the center of the initial
% image mask
stats = regionprops(MaskInit,'centroid');
centroid = stats.Centroid;
for i = 1:stacklength
    A_mask(:,:,i) = A(:,:,i);
    A_mask(:,:,i) = logical(MaskInit);
    A_region(:,:,i) = A_bilatfilt(:,:,i).*A_mask(:,:,i);
    A_region = A_region(:,:,i);
    [Max,MaxInd] = max(A_region,[],'all','linear');
    [MaxIndx,MaxIndy] = ind2sub([size(A,1) size(A,2)],MaxInd);
    distances = sqrt(sum(bsxfun(@minus, [MaxIndx,MaxIndy], centroid).^2,2));              % compute the Euclidean distances
    closest = MaxInd(find(distances==min(distances)),:);
    [x(i),y(i)] = ind2sub([size(A,1) size(A,2)],MaxInd);
    MaxVal(i) = Max(find(distances==min(distances)),:);
end
MaxVal = MaxVal';
x = x';
y = y';


%% Now exclude seed points that are below a threshold value
for i = 1:stacklength
    if MaxVal(i) >= 0
        T(i) = 1;
    else 
        T(i) = 0;
    end
end
T = T';

for i = 1:stacklength
    x(i) = x(i)*T(i);
    y(i) = y(i)*T(i);
end


%% Conduct seed growing algorithm on all image slices excluding ones
%  without seed points
for i = 1:stacklength
    if x(i) == 0 & y(i) == 0
        figure
        imshow(A(:,:,i))
        F(i) = getframe;
        A_final(:,:,:,i)=F.cdata;
        A_finalS(:,:,:,i) = frame2im(F(i));
        close
        continue
    end
        [final, colored, pouthisteq] = Image_Analysis_Fxn(A(:,:,i),x(i),y(i),tol);   % run image analysis fxn and store final images
        disc(:,:,i) = final(:,:);                       % store analysis for each slice
        J(:,:,:,i) = colored(:,:,:);                    % store labeloverlay for each slice
        figure
        imshow(colored)
        F(i) = getframe;
        A_final(:,:,:,i)=F.cdata;
        A_finalS(:,:,:,i) = frame2im(F(i));
        ImgPre(:,:,i) = pouthisteq(:,:);
        close
end

% figure
% montage(ImgPre)

% figure
% montage(A_finalS)                                          % display all analyzed image slices from stack
% title('Analysis Overlayed with Preprocessed Image')


%% If segmented image touches border disregard all parts of segmented 
% image boundaries that fall outside of the initial image mask dilated
MaskInit = imclose(MaskInit,strel('disk',8));

for i = 1:stacklength
    if x(i) == 0 & y(i) == 0
        figure
        imshow(A(:,:,i))
        F(i) = getframe;
        A_final(:,:,:,i)=F.cdata;
        A_finalS(:,:,:,i) = frame2im(F(i));
        close
    continue
    end
    ClearBord(:,:,i) = imclearborder(disc(:,:,i));             % if disc touches border clear
    true = find(ClearBord(:,:,i));                      % find images that have been cleared
    if isempty(true) == 1
        disc(:,:,i) = disc(:,:,i).*logical(MaskInit);
    end
end


%% If segmented image touches image boundary then exclude that disc
for i = 1:stacklength
    if x(i) == 0 & y(i) == 0
        figure
        imshow(A(:,:,i))
        F(i) = getframe;
        A_final(:,:,:,i)=F.cdata;
        A_finalS(:,:,:,i) = frame2im(F(i));
        close
        continue
    end
    ClearBord(:,:,i) = imclearborder(disc(:,:,i));             % if disc touches border clear
    true = find(ClearBord(:,:,i));                      % find images that have been cleared
    if isempty(true) == 1
        disc(:,:,i) = zeros(size(disc(:,:,i)));
        x(i) = 0;
        y(i) = 0;
    end 
end


%% Disregard segmented images that are not within disc boundaries obtained
% from GUI
for i = 1:stacklength
    if i<discbegin || i>discend
        disc(:,:,i) = zeros(size(A(:,:,i)));
    end
end


%% Use imclose to fill holes in the disc segments
for i = 1:stacklength
    disc(:,:,i) = imclose(disc(:,:,i),strel('disk',5));
end


%% Compare analyzed images to original images without processing
for i = 1:stacklength
    if x(i) == 0 & y(i) == 0
        figure;
        imshow(A(:,:,i));
        F(i) = getframe;
        ImgComp(:,:,:,i)=F.cdata;
        ImgComp(:,:,:,i) = frame2im(F(i));
        close
        continue
    end
    if i <discbegin || i>discend
               figure;
        imshow(A(:,:,i));
        F(i) = getframe;
        ImgComp(:,:,:,i)=F.cdata;
        ImgComp(:,:,:,i) = frame2im(F(i));
        close
        continue 
    end
    ImgComp(:,:,:,i) = labeloverlay(A(:,:,i),disc(:,:,i));
    figure;
    imshow(ImgComp(:,:,:,i));
    hold on
    plot(x(i),y(i),'+','MarkerFaceColor','red')
    set(gca, 'YDir', 'reverse');
    
    F(i) = getframe;
    ImgComp(:,:,:,i)=F.cdata;
    ImgComp(:,:,:,i) = frame2im(F(i));
    close
end
% figure
% montage(ImgComp)
% title('Seed Analysis Overlayed with Original Raw Image')


%% Export segmented images into a .tif file

% Here we reformat the isolated disc into uint8 image format
for i = 1:stacklength
    discnew(:,:,i)=im2uint8(disc(:,:,i));
end

finalfile = append('FinalDisc','_',ImgStack);
imwrite(discnew(:,:,1),finalfile);
for i = 2:stacklength
    imwrite(discnew(:,:,i),finalfile,'WriteMode','append','Compression','none');
end



final_disc = discnew;
final_overlay = ImgComp;

end