% Sample use of PointCloud2Image(...)
% 
% The following variables are contained in the provided data file:
%       BackgroundPointCloudRGB,ForegroundPointCloudRGB,K,crop_region,filter_size
% None of these variables needs to be modified


clc
clear all
% load variables: BackgroundPointCloudRGB,ForegroundPointCloudRGB,K,crop_region,filter_size)
load data.mat

data3DC = {BackgroundPointCloudRGB,ForegroundPointCloudRGB};
R       = eye(3);

move    = [0 0 -0.025]'; % Choosing a value which does not go beyond the foreground

for step = 0:74
   tic
   fname       = sprintf('SampleOutput3%03d.jpg',step);
   display(sprintf('\nGenerating %s',fname));
   t           = step * move;
   Z    =  3.8650 + t(3); % Gives a bounding box for the foreground as approximately 640x400 (actual value -> 627 x 438)...
                          % 3.8650 is found out approximately by finding
                          % the value of t which approximately bounds the
                          % foreground to the bounding box
   ratio_x = 2759.5/3.8650;
   ratio_y = 2764.2/3.8650;
   K(1,1) = ratio_y * Z;
   K(2,2) = ratio_x * Z;
%    K(1,1)          = (640/(0.0529 + 0.5)) * Z; % fy
 
%    K(2,2)       = (400/(0.3287 + 0.0589)) * Z;  % fx
   M           = K*[R t];
   im          = PointCloud2Image1(M,data3DC,crop_region,filter_size);
   imwrite(im,fname);
   toc
end