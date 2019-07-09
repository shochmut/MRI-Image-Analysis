% Skyler Hochmuth
% Colorado State University
% Walter Scott School of Biomedical Engineering
% Spring 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab file converts the pixelated MRI images into a coordinate    %
% system. Images are then analyzed using a peak analysis method to find   % 
% the pixels that contain the nucleus pulposis (white) contained within   %
% the annulus fibrosis (black)                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% First add the subpath of the directory: (You must be at the root of the
% directory)
addpath(genpath('.'));


%% Now we load the image file (image stack containing multiple slice images)
% that is in tif format
% ImgStack = '19_Standing_L2L3.tif';                  % set up file name
% ImgStack = '14_Standing_L5S.tif';
ImgStack = '14_Standing_L2L3.tif';
t = Tiff(ImgStack,'r');                             % open and read tif file


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
for i = 13
    A_bilatfilt = imbilatfilt(A(:,:,i));
    A_histeq = histeq(A(:,:,i));
    sharp = imsharpen(A(:,:,i));
    hist_bilatfilt = histeq(A_bilatfilt);
    sharp_hist_bilatfilt = imsharpen(hist_bilatfilt,'Amount',1,'Radius',2.5);
end
figure
montage({A(:,:,i),A_bilatfilt,A_histeq,sharp_hist_bilatfilt},'Size',[1 4])   % display the 3 images    


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
tol = 0.2;

%% Run the seed growing algorithm on initial image
[final, colored, pouthisteq] = Image_Analysis_Fxn(InitImg,x,y,tol); % run image analysis fxn and store final images
    MaskInit(:,:) = final(:,:);                       % store inital image mask
figure
imshow(colored)
title('Initial Mask Generation')


%% Use the obtained initial mask to find seed points in all other image slices
for i = 1:stacklength
    A_mask(:,:,i) = A(:,:,i);
    A_mask(:,:,i) = logical(MaskInit);
    A_region(:,:,i) = A(:,:,i).*A_mask(:,:,i);
    A_region = A_region(:,:,i);
    [Max,MaxInd] = max(A_region(:));
    [x(i),y(i)] = ind2sub([size(A,1) size(A,2)],MaxInd);
    MaxVal(i) = Max;
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

figure
montage(ImgPre)

figure
montage(A_finalS)                                          % display all analyzed image slices from stack
title('Analysis Overlayed with Preprocessed Image')


%% If segmented image touches image boundary then imopen that disc
% 
% for i = 1:stacklength
%     if x(i) == 0 & y(i) == 0
%         figure
%         imshow(A(:,:,i))
%         F(i) = getframe;
%         A_final(:,:,:,i)=F.cdata;
%         A_finalS(:,:,:,i) = frame2im(F(i));
%         close
%         continue
%     end
%     ClearBord(:,:,i) = imclearborder(disc(:,:,i));             % if disc touches border clear
%     true = find(ClearBord(:,:,i));                      % find images that have been cleared
%     if isempty(true) == 1
%         disc(:,:,i) = imopen(disc(:,:,i),strel('disk',3));
%         disc(:,:,i) = imclearborder(disc(:,:,i));
%     end 
% end

%% If segmented image touches image boundary then lower tolerance on that
% disc region growing algorithm
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
%     if tol <= 0
%         x(i) = 0; 
%         y(i) = 0; 
%         continue
%     end
    ClearBord(:,:,i) = imclearborder(disc(:,:,i));             % if disc touches border clear
    true = find(ClearBord(:,:,i));                      % find images that have been cleared
    if isempty(true) == 1
        tol = 0.4
        for j = 20:1
            tol = tol - (tol/j);
            [final, colored, pouthisteq] = Image_Analysis_Fxn(A(:,:,i),x(i),y(i),tol);   % run image analysis fxn and store final images
            disc(:,:,i) = final;                       % store analysis for each slice
            J(:,:,:,i) = colored;                    % store labeloverlay for each slice
        end 
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


%% Compare analyzed images to original images without processing
for i = 1:stacklength;
    if x(i) == 0 & y(i) == 0
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

    F(i) = getframe;
    ImgComp(:,:,:,i)=F.cdata;
    ImgComp(:,:,:,i) = frame2im(F(i));
    close
end
figure
montage(ImgComp)
title('Seed Analysis Overlayed with Original Raw Image')


%% Export segmented original images into tif file to view on fiji
imwrite(ImgComp(:,:,:,1),'19_Standing_L2L3_Segmented.tif')
for i = 2:stacklength;
    imwrite(ImgComp(:,:,:,i),'19_Standing_L2L3_Segmented.tif',...
        'WriteMode','append')
end

for i = 1:stacklength
    discnew(:,:,i)=im2uint8(disc(:,:,i));
end
imwrite(discnew(:,:,1),'19_Standing_L2L3_Isolated.tif')
for i = 2:stacklength;
    imwrite(discnew(:,:,i),'19_Standing_L2L3_Isolated.tif',...
        'WriteMode','append')
end







% %% Try canny edge detection method as an alternative to seed-grow algorithm
% 
% % Pre process the image
% for i = 1:stacklength
%     A_canny = A(:,:,i);
%     A_process = imbilatfilt(A_canny);                         % apply an edge-preserving bilateral filter
%     pouthisteq = histeq(A_canny);                            % apply histogram eq to enhance contrast
%     sharp = imsharpen(pouthisteq);                       % sharpen the filtered image
%     pouthisteq = sharp;
%     
% %     A_canny = sharp;                                          % set the image equal to the sharpened image
%     pouthisteq = histeq(A_canny);
% 
%     figure
%     BW2 = edge(pouthisteq,'Canny',0.5,1);            % conduct canny edge processing
%     L_total = logical(BW2);                             % convert processed image to logical
%     Segout = A_canny;
%     Segout(L_total) = 255;                              % define edge boundaries
%     imshow(Segout)                                      % display canny edge 
%     title('Canny Edge Detector')
% 
%     BW2 = imfill(BW2,4,'holes');
%     
%     L_total = logical(BW2);
% 
%     se = strel('disk',2);
%     BW2 = imopen(BW2,se);
%     BW2 = bwareaopen(BW2,50);
%     BW2 = imclose(BW2,strel('disk',3));
%     
%     figure
%     imshow(labeloverlay(A(:,:,i),BW2))
%     F(i) = getframe;
%     ImgCanny(:,:,:,i)=F.cdata;
%     ImgCanny(:,:,:,i) = frame2im(F(i));
%     close
% end
% 
% figure
% montage(ImgCanny)
% title('Canny Image Processing')
% 
% 










% %% Use a GUI to determine the seed point for the image slices
% 
% figure
% x = []; y = [];                                     % preallocate empty points
% for i = 1:stacklength                               % initiate for loop
%     if isempty(x) == 1 && isempty(y) == 1           % if parameters are empty
%         imshow(imread(ImgStack,i));
%         text(0,-10,'Double click on seed point')
%         text(0,-5,'If no visible seed point press enter')
%         [x,y] = getpts;                             % get seed point
%     end
%     if isempty(x) == 0 && isempty(y) == 0           % if seed points chosen
%         break                                       % break for loop
%     end
% end
% leftslice = i;
% x = round(x);                                       % round seed point
% y = round(y);
% close all                                           % close images
% 
% 
% 
% %% Now run the image analysis function for each image slice
% 
% % Run a for loop to iterate through each slice image
% for i = 1:stacklength;
%     A = imread(ImgStack,i);                         % read image slice
%     [final, colored] = Image_Analysis_Fxn(A,x,y);   % run image analysis fxn and store final images
%     disc(:,:,i) = final(:,:);                       % store analysis for each slice
%     J(:,:,:,i) = colored(:,:,:);                    % store labeloverlay for each slice
%     F(i) = getframe;
%     A_final(:,:,:,i)=F.cdata;
%     A_final(:,:,:,i) = frame2im(F(i));
% end
% figure
% montage(A_final)                                          % display all analyzed image slices from stack
% title('Analysis Overlayed with Preprocessed Image')
% 
% %% Compare analyzed images to original images without processing
% %for i = 1:stacklength;
% for i = 1:stacklength;
%     A = imread(ImgStack,i);
%     ImgComp(:,:,:,i) = labeloverlay(A,disc(:,:,i));
% end
% figure
% montage(ImgComp)
% title('Analysis Overlayed with Original Raw Image')
%     
%     
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
