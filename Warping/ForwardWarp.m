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

p1_prime = [-x1 -y1 -1 0 0 0 x1*x11 y1*x11 x11; 0 0 0 -x1 -y1 -1 x1*y11 y1*y11 y11];
p2_prime = [-x2 -y2 -1 0 0 0 x2*x22 y2*x22 x22; 0 0 0 -x2 -y2 -1 x2*y22 y2*y22 y22];
p3_prime = [-x3 -y3 -1 0 0 0 x3*x33 y3*x33 x33; 0 0 0 -x3 -y3 -1 x3*y33 y3*y33 y33];
p4_prime = [-x4 -y4 -1 0 0 0 x4*x44 y4*x44 x44; 0 0 0 -x4 -y4 -1 x4*y44 y4*y44 y44];

P = [p1_prime; p2_prime; p3_prime; p4_prime];

[U,S,V] = svd(P);
% 
H = V(:,end); %Homography Matrix

H = vec2mat(H,3);
% H_inv = inverse(H);

% H = [3 0; 0 3];
% H(3,3) = 1;

X = cell(1,2);
X_new = zeros(3,1);
% 
k = 1;
for i = 1:rows
    for j = 1:cols
        X{1,k} =  H*[j;i;1];
        k = k+1;
    end
end
H_inverse = inv(H);
H_inv = H_inverse;
% disp(k)

X_prime = cell(1,2);
x_X_prime_skipped = zeros(2,1);
y_X_prime_skipped = zeros(2,1);

x_X_prime = zeros(2,1);
y_X_prime = zeros(2,1);

%% Transformed co-ordinates

kk = 1;

for kk = 1:k-1
    if round(X{1,kk}(1,1)/X{1,kk}(3,1)) <= 0 || round(X{1,kk}(2,1)/X{1,kk}(3,1))<= 0 || round(X{1,kk}(1,1)/X{1,kk}(3,1)) > 940 || round(X{1,kk}(2,1)/X{1,kk}(3,1)) > 500
%                 I_warped(round(abs(X{1,kk}(1,1)/X{1,kk}(3,1))), abs(round(X{1,kk}(2,1)/X{1,kk}(3,1))),:) = 0;
        continue;
    else
        X_prime{1,kk}(1,1) = round(X{1,kk}(1,1)/X{1,kk}(3,1));
        X_prime{1,kk}(2,1) = round(X{1,kk}(2,1)/X{1,kk}(3,1));
        
        x_X_prime(kk,1) = round(X{1,kk}(1,1)/X{1,kk}(3,1));
        y_X_prime(kk,1) = round(X{1,kk}(2,1)/X{1,kk}(3,1));
       
%         I_warped(round(abs(X{1,kk}(1,1)/X{1,kk}(3,1))),round(abs(X{1,kk}(2,1)/X{1,kk}(3,1))), :) = I(i,j,:);
%             kk = kk+1;    
    end
%     end
end

x_X_prime = x_X_prime(x_X_prime ~= 0);
y_X_prime = y_X_prime(y_X_prime ~= 0);
XX_prime = [x_X_prime, y_X_prime];



% X_prime =  X_prime(~cellfun('isempty',X_prime));


c = 0;
ll = 1;
for i = 1:rows
    for j = 1:cols
        if isempty(X_prime{1,ll})
           ll = ll + 1;
           continue;
        else
           I_warped(X_prime{1,ll}(1,1),X_prime{1,ll}(2,1),:) = I(i,j,:);
           ll = ll+1;
        end            
        if ll > length(X_prime)
           c = c+1;
           break;
        end      
    end
        if c == 1
           break;
        end
end
%



Red = I_warped(:,:,1) == 0;
Green = I_warped(:,:,2) == 0;
Blue = I_warped(:,:,3) == 0;

RGB_black = Red.*Green.*Blue; 

[x_skipped, y_skipped] = find(RGB_black == 1);
X_skipped = [x_skipped, y_skipped];

% X = zeros(2,1);
% [X] = list(X);

% X_skipped_x = zeros(1,2);
% X_skipped_y = zeros(1,2);
% rr = 1;
% for i = 1:length(X_prime)
%     for j = 1:length(X)
%         if isequal(X(j,:), X_prime(i,:))
%            continue;
%         else
%             X_skipped_x(rr) = i;
%             X_skipped_y(rr) = j;
%             rr = rr+1;
%         end
%     end
% end



% 
% ll = 1;
% o = 1;
% X_skipped = zeros(1,2);
% for i = 1:size(I_warped,1)
%     for j = 1:size(I_warped,2)
%         Y = [i,j];
%         if find (XX_prime(:,:) == Y) >= 1
%            continue;
%         else
%             X_skipped(ll,1) = i;
%             X_skipped(ll,2) = j;
%             ll = ll+1;
%         end
%     end
% end



X_new_cord = cell(1,2);
for i = 1:length(X_skipped)
        X_new =  H_inv*[X_skipped(i,1);X_skipped(i,2);1];
        X_new_cord = [X_new(1,1)/X_new(3,1), X_new(2,1)/X_new(3,1)];
        R11_red = I( floor(X_new_cord(2)),floor(X_new_cord(1)),1);
        R12_red = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),1);
        R22_red = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,1);
        R21_red = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,1);
        
        R11_green = I( floor(X_new_cord(2)),floor(X_new_cord(1)),2);
        R12_green = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),2);
        R22_green = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,2);
        R21_green = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,2);
        
        R11_blue = I( floor(X_new_cord(2)),floor(X_new_cord(1)),3);
        R12_blue = I( floor(X_new_cord(2))+1,floor(X_new_cord(1)),3);
        R22_blue = I( floor(X_new_cord(2))+1,floor(X_new_cord(1))+1,3);
        R21_blue = I( floor(X_new_cord(2)),floor(X_new_cord(1))+1,3);
        R = (X_new_cord - floor(X_new_cord));
        a = R(1); b = R(2);
        
        I_warped(X_skipped(i,1),X_skipped(i,2),1) = round((1-a)*(1-b)*R11_red + a*(1-b)*R21_red + a*b*R22_red + b*(1-a)*R12_red);
        I_warped(X_skipped(i,1),X_skipped(i,2),2) = round((1-a)*(1-b)*R11_green + a*(1-b)*R21_green + a*b*R22_green + b*(1-a)*R12_green);
        I_warped(X_skipped(i,1),X_skipped(i,2),3) = round((1-a)*(1-b)*R11_blue + a*(1-b)*R21_blue + a*b*R22_blue + b*(1-a)*R12_blue);
%         e = e+1;
end





I_warped = uint8(I_warped);
figure
imshow(I_warped);

 