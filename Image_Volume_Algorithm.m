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

% First add the subpath of the directory: (You must be at the root of the
% directory)
addpath(genpath('.'));


%% Now we load the image and crop it

A = imread('2_Standing_L2L311.tif');                % read image file
figure                                              % create figure
imshow(A)                                           % show original image
axis on
A_crop = A(20:45,10:55);                            % crop image
figure                                              % create another figure    
imshow(A_crop)                                      % show cropped image

x=1:length(A_crop);                                % narrow x-axis to cropped length

zero = zeros(size(A),'like',A);                                             %Pre-allocate matrix of zeros the same size as image data
zero1 = zero;
zero2 = zero;

figure
title('Horizonal Peak Analysis')
ylabel('Image Grayscale Value')
xlabel('X Pixel Coordinate')
for y = 1:length(A_crop(:,1))                      % y axis, narrows image down   
    hold on
      findpeaks(double(A_crop(y,:)),'Annotate','peaks','MinPeakHeight',...  %find all peaks above threshold and plot peak analysis
      100);                                                                  %min peak height chosen to exclude noise
      ylim([0 250])

      [pks,locs,~,~] = findpeaks(double(A_crop(y,:)),'Annotate',...         %create vector of peak heights and locations
          'peaks','MinPeakHeight', 100);
      localmin = islocalmin(A_crop(y,:));                                   %find all local min's indexed positions            
      localmin_x = x(localmin);                                             %create vector of local min x locations
      
      if isempty(locs)                                                      %if no peaks are found force next iteration
          right(y) = 0;
          left(y) = 0;
          continue
      end
      
      leftpeak_x = locs(1);                                                 %find farthest peak to left
      rightpeak_x = locs(end);                                              %find farthest peak to right
      
      [~,leftbound_I] = find(localmin_x<leftpeak_x,1,'last');               %find index position of local min closest and to left of left peak
      [~,rightbound_I] = find(localmin_x>rightpeak_x,1,'first');            %find index position of local min closest and to right of right peak
      
      leftbound = localmin_x(leftbound_I);                                  %determine leftbound x position
      rightbound = localmin_x(rightbound_I);                                %determine rightbound x position
      
      if isempty(leftbound)                             %%%%%%%%%%%%%%%%%
          continue                                %%%%%%to make more accurate, change these to use next innermost peak
      end                                         %%%%%%%to find bounds
                                                        %%%%%%%%%%%%%%%%%
      if isempty(rightbound)
          continue
      end
      
      plot(x(leftbound),double(A_crop(y,leftbound)),'g*')                   %plot bounds 
      plot(x(rightbound),double(A_crop(y,rightbound)),'g*')
      right(y) = rightbound;                                                %store rightbound value
      left(y) = leftbound;                                                  %store leftbound value
    
      width(y) = rightbound - leftbound;                                    %sum width
      zero1(20+y,(right(y)+10)) = A(20+y,(right(y)+10));                               %sub A values for right bounds 
      zero1(20+y,(left(y)+10)) = A(20+y,(left(y)+10));                                 %sub A values for left bounds
end

Disc = sum(width)                                                           %sum all widths to get disc area in [pixels]

figure                                                                 
L = logical(zero1);
BWoutline1 = bwperim(L);
Segout1 = A; 
Segout1(BWoutline1) = 255; 
imshow(Segout1)
title('Horizonal Peak Analysis')


y=1:length(A_crop(:,1));                                % x axis narrowed


figure
title('Vertical Peak Analysis')
ylabel('Image Grayscale Value')
xlabel('Y Pixel Coordinate')
for x = 1:length(A_crop(1,:))                      % x axis, narrows image down   
    hold on
      findpeaks(double(A_crop(:,x)),'Annotate','peaks','MinPeakHeight',...  %find all peaks above threshold and plot peak analysis
      100);                                                                  %min peak height chosen to exclude noise
      ylim([0 250])

      [pks,locs,~,~] = findpeaks(double(A_crop(:,x)),'Annotate',...         %create vector of peak heights and locations
          'peaks','MinPeakHeight', 100);
      localmin = islocalmin(A_crop(:,x));                                   %find all local min's indexed positions            
      localmin_y = y(localmin);                                             %create vector of local min y locations
      
      if isempty(locs)                                                      %if no peaks are found force next iteration
          top(x) = 0;
          bottom(x) = 0;
          continue
      end
      
      bottompeak_y = locs(1);                                                 %find farthest peak to bottom
      toppeak_y = locs(end);                                              %find farthest peak to top
      
      [~,bottombound_I] = find(localmin_y<bottompeak_y,1,'last');               %find index position of local min closest and to left of left peak
      [~,topbound_I] = find(localmin_y>toppeak_y,1,'first');            %find index position of local min closest and to right of right peak
      
      bottombound = localmin_y(bottombound_I);                                  %determine bottombound y position
      topbound = localmin_y(topbound_I);                                %determine topbound y position
      
      if isempty(bottombound)                             %%%%%%%%%%%%%%%%%
          continue                                %%%%%%to make more accurate, change these to use next innermost peak
      end                                         %%%%%%%to find bounds
                                                        %%%%%%%%%%%%%%%%%
      if isempty(topbound)
          continue
      end
      
      plot(y(bottombound),double(A_crop(bottombound,x)),'g*')                   %plot bounds 
      plot(y(topbound),double(A_crop(topbound,x)),'g*')
      top(x) = topbound;                                                %store rightbound value
      bottom(x) = bottombound;                                                  %store leftbound value
    
      height(x) = topbound - bottombound;                                    %sum width
      zero2((top(x)+20),x+10) = A((top(x)+20),x+10);                               %sub A values for right bounds 
      zero2((bottom(x)+20),x+10) = A((bottom(x)+20),x+10);                                 %sub A values for left bounds
end


figure                                                                 
L = logical(zero2);
BWoutline2 = bwperim(L);
Segout2 = A; 
Segout2(BWoutline2) = 255; 
imshow(Segout2)
title('Vertical Peak Analysis')


figure
L1 = logical(zero1);
L2 = logical(zero2);
L = zero1+zero2;
L_total = logical(L);
Segout = A;
Segout(L_total) = 255;
imshow(Segout)
title('Vertical and Horizontal Analysis Combined')


index = 0;
for i = 1:65
    for j = 1:65
        if L_total(i,j) == 1
            index = index + 1;
            X(index) = j;
            Y(index) = i;
        end
    end
end

% % Now drop any rows with zeros 
Tf1 = X==0;
Tf2 = Y==0;
X(Tf1)=[];
Y(Tf2)=[];
XY = [X',Y'];


figure
scatter(X,Y)
set(gca,'Ydir','reverse')


figure
Ellipse = fit_ellipse(X,Y);                                                 % Run ellipse fitting function


figure
hold on
scatter(X,Y)
set(gca,'Ydir','reverse')
ellipse(Ellipse.a,Ellipse.b,-Ellipse.phi,Ellipse.X0_in,Ellipse.Y0_in)

% Now we plot the final formed ellipse onto the MRI disc image to see the
% accuracy of our results.

figure
imshow(Segout)
hold on
axis on
axis tight
ellipse(Ellipse.a,Ellipse.b,-Ellipse.phi,(Ellipse.X0_in),Ellipse.Y0_in) 




% Now to conduct spline interpolation to see if we can get a better fit to
% the nucleus pulposis compared to the ellipse fit method
% First find the point furthest to the left and point furthest to the right
[left_spline,left_index] = min(XY(:,1));                                    %find left point and index
[right_spline,right_index] = max(XY(:,1));                                  %find right point and index

% Now to isolate data above these points for top spline curve
XY_topspline = zeros(size(XY));                                             %Pre-allocate matrix for speed
for i = 1:length(Y)
    if XY(i,2) <= XY(left_index,2)                                           %Select for data points above left point
        XY_topspline(i,:) = XY(i,:);
    end
end
XY_topspline = [nonzeros(XY_topspline(:,1)),nonzeros(XY_topspline(:,2))];   %Get rid of preallocated zeros
XY_topspline = sortrows(XY_topspline,1);                                   %Rearrange data points in X ascending order
pt = interparc(100,XY_topspline(:,1),XY_topspline(:,2),'spline');           %Distance based interpolation along general curve in space function

%Now plot the interpolation compared with data points
figure
plot(XY_topspline(:,1),XY_topspline(:,2),'r*',pt(:,1),pt(:,2),'b-o')
        
%Now plot the interpolation onto MRI image
figure
imshow(Segout)
hold on
axis on
axis tight
plot(pt(:,1),pt(:,2),'b','linewidth',3)

% Now to isolate the points under the leftmost point for second spline
% curve
XY_bottomspline = zeros(size(XY));                                             %Pre-allocate matrix for speed
for i = 1:length(Y)
    if XY(i,2) > XY(left_index,2)                                           %Select for data points above left point
        XY_bottomspline(i,:) = XY(i,:);
    end
end
XY_bottomspline = [nonzeros(XY_bottomspline(:,1)),nonzeros(XY_bottomspline(:,2))];   %Get rid of preallocated zeros
XY_bottomspline = [XY(left_index,:); XY_bottomspline];
XY_bottomspline = sortrows(XY_bottomspline,1);                                   %Rearrange data points in X ascending order
spline_bottom = interparc(100,XY_bottomspline(:,1),XY_bottomspline(:,2),'spline');           %Distance based interpolation along general curve in space function

%Now plot the interpolation compared with data points
figure
plot(XY_bottomspline(:,1),XY_bottomspline(:,2),'r*',spline_bottom(:,1),spline_bottom(:,2),'b-o')

%Now plot the interpolation onto MRI image
figure
imshow(Segout)
hold on
axis on
axis tight
plot(spline_bottom(:,1),spline_bottom(:,2),'b','linewidth',3)

%Now put bottom and top interpolations together and onto MRI image
figure
imshow(Segout)
hold on
axis on
axis tight
plot(spline_bottom(:,1),spline_bottom(:,2),'g','linewidth',2)
plot(pt(:,1),pt(:,2),'g','linewidth',2)


%% Try spline arc interpolation for periodic closed curve instead
XY_topspline = sortrows(XY_topspline,1,'descend');
XY_arc = [XY_bottomspline;XY_topspline];
spline_arc = interparc(200,XY_arc(:,1),XY_arc(:,2),'csape');
figure
plot(XY_arc(:,1),XY_arc(:,2),'r*',spline_arc(:,1),spline_arc(:,2),'b-o')
figure
imshow(Segout)
hold on
axis on
axis tight
plot(spline_arc(:,1),spline_arc(:,2),'g','linewidth',2)

% Sort the data points into clockwise arrangment










% % Now to interpolate between the points to form a final disc (ellipse)
% % First we form a matrix containing the X (column 1) and Y values (column 2)
% % for the points obtained from the peak analysis
% XY1 = [(1:46)',top'];                                                       % Convert left, right, top,
% XY2 = [(1:46)',bottom'];                                                    % and bottom coordinate values
% XY3 = [left',(1:26)'];                                                      % into XY coordinates 
% XY4 = [right',(1:26)'];
% 
% % Now drop any rows with zeros 
% Tf1 = XY1(:,2)==0;
% Tf2 = XY2(:,2)==0;
% Tf3 = XY3(:,1)==0;
% Tf4 = XY4(:,1)==0;
% XY1(Tf1,:)=[];
% XY2(Tf2,:)=[];
% XY3(Tf3,:)=[];
% XY4(Tf4,:)=[];
% 
% % Now concatenate the matrices together to get all points in one matrix
% XY = [XY1;XY2;XY3;XY4];
% X = XY(:,1);                                                                % Separate X points into X vector
% Y = XY(:,2);                                                                % Separate Y points into Y vector

% figure
% Ellipse = fit_ellipse(X,Y);                                                 % Run ellipse fitting function
% h = ellipse(Ellipse.a,Ellipse.b,Ellipse.phi,Ellipse.X0,Ellipse.Y0)          % Plot ellipse using ellipse function
% axis([0 46 0 26])                                                           % Set axis of plot


% Now we plot the final formed ellipse onto the MRI disc image to see the
% accuracy of our results.

% figure
% imshow(Segout)
% hold on
% axis on
% axis tight
% ellipse(Ellipse.a,Ellipse.b,Ellipse.phi,(Ellipse.X0),Ellipse.Y0) 





% K = convhull(X,Y);                                                          % Take the convex hull
% plot(X(K),Y(K),'g')




% figure
% zero1 = zeros(65);
% L1 = zero1(X,Y)+100;
% % L1 = logical(zero1);
% CH = bwconvhull(L1);
% Segout = A;
% Segout(L_total) = 255;
% imshow(CH)
% title('Convex Hull of Data Points')











