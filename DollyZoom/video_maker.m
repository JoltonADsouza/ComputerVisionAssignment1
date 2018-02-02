clc
close all;
clear all;
 
images = cell(1,75); % Create a cell of images

   srcFiles = dir('C:\Users\Jolton\Desktop\Images\*.jpg'); % Read the directory of the created images
 for i = 1: length(srcFiles)
  filename = strcat('C:\Users\Jolton\Desktop\Images\',srcFiles(i).name);
  images{i} = imread(filename);
 end

 writerObj = VideoWriter('Output_video.wmv');
 writerObj.FrameRate = 15; % 75 Images, 15 images per second gives a 5 second video for 75 images

 open(writerObj); % Open video writer
 
 for u=1:length(images)
     frame = im2frame(images{u});
     writeVideo(writerObj, frame);
 end
 close(writerObj);  % close the writer object