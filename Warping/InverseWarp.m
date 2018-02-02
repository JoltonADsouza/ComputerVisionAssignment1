clc
clear all;
close all;

I = imread('C:\Users\Jolton\Desktop\Files\Homeworks\Semester3\3D Computer Vision\Homework1\basketball-court.ppm');
I_warped = zeros(940,500,3);
% I_warped(:,:,1) = 255; I_warped(:,:,2) = 255; I_warped(:,:,3) = 0;

[rows,cols,color] = size(I);
figure
imshow(I);

%% Finding the Homography matrix

x1 = 245; y1 = 45;
x2 = 420; y2 = 70;
x3 = 1; y3 = 195;
x4 = 300; y4 = 295;

x11 = 1; y11 = 1;
x22 = 1; y22 = 500;
x33 = 940; y33 = 1;
x44 = 940; y44 = 500;
% 
% x1 = 0; y1 = 194;
% x2 = 245 ; y2 = 45;
% x3 = 418; y3 = 70;
% x4 = 296; y4 = 300;

% x11 = 0; y11 = 0;
% x22 = 939 ; y22 = 0;
% x33 = 939; y33 = 499;
% x44 = 0; y44 = 499;

p1_prime = [-x1 -y1 -1 0 0 0 x1*x11 y1*x11 x11; 0 0 0 -x1 -y1 -1 x1*y11 y1*y11 y11];
p2_prime = [-x2 -y2 -1 0 0 0 x2*x22 y2*x22 x22; 0 0 0 -x2 -y2 -1 x2*y22 y2*y22 y22];
p3_prime = [-x3 -y3 -1 0 0 0 x3*x33 y3*x33 x33; 0 0 0 -x3 -y3 -1 x3*y33 y3*y33 y33];
p4_prime = [-x4 -y4 -1 0 0 0 x4*x44 y4*x44 x44; 0 0 0 -x4 -y4 -1 x4*y44 y4*y44 y44];

P = [p1_prime; p2_prime; p3_prime; p4_prime];

[U,S,V] = svd(P);
% 
H = V(:,end); %Homography Matrix

H = vec2mat(H,3);

H_inv = inv(H);

X = zeros(3,1);
X_new_cord = zeros(1,2);

a = 0;
b = 0;

% 
for i = 1:size(I_warped,2)
    for j = 1:size(I_warped,1)
% for i = 1
%     for j = 2
            X =  H_inv*[j;i;1];
        
           
            
            X_new_cord = [X(1,1)/X(3,1), X(2,1)/X(3,1)];
            
            % Bilinear Interpolation
            R = (X_new_cord - floor(X_new_cord));
                if uint8(R(1)) == 1 && uint8(R(2)) == 1
%                      a = 0; b = 0;
                     
                     I_warped(j,i,1) = I(uint32(X_new_cord(2)),uint32(X_new_cord(1)),1);
                     I_warped(j,i,2) = I(uint32(X_new_cord(2)),uint32(X_new_cord(1)),2);
                     I_warped(j,i,3) = I(uint32(X_new_cord(2)),uint32(X_new_cord(1)),3);
                     
                else
                    a = R(1); b = R(2);
                    if floor(X_new_cord(1)) < 1 || floor(X_new_cord(1)) > cols || floor(X_new_cord(2)) < 1 || floor(X_new_cord(2)) > rows
                        continue;
                    else
                     R11 = I( floor(X_new_cord(2)),floor(X_new_cord(1)),1); R12 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),1);
                     R22 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,1); R21 = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,1);
        
                     G11 = I( floor(X_new_cord(2)),floor(X_new_cord(1)),2); G12 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),2);
                     G22 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,2); G21 = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,2);
        
                     B11 = I( floor(X_new_cord(2)),floor(X_new_cord(1)),3); B12 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),3);
                     B22 = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,3); B21 = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,3);
                    
                     I_warped(j,i,1) = round((1-a)*(1-b)*R11 + a*(1-b)*R21 + a*b*R22 + b*(1-a)*R12);
                     I_warped(j,i,2) = round((1-a)*(1-b)*G11 + a*(1-b)*G21 + a*b*G22 + b*(1-a)*G12);
                     I_warped(j,i,3) = round((1-a)*(1-b)*B11 + a*(1-b)*B21 + a*b*B22 + b*(1-a)*B12);
                    end
                end  
    end     
end
I_warped = uint8(I_warped);
I_warped = imrotate(I_warped, -90);
figure
imshow(I_warped);