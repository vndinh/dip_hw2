clear; clc;

Problem2();

function Problem2()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw2\Test_images');
file = fopen(fullfile(imgdir,'\Color_monarch_512x512.raw'),'rb');
color_monarch = fread(file,fliplr([512,512*3]),'*uint8')';
fclose(file);

%%
%%---------------------Insert code below ----------------------%%
[R, G, B, rgb_monarch] = changeRGB(color_monarch);
[H, S, I, hsi_monarch] = RGB2HSI(rgb_monarch);
[Y, Cb, Cr, ycbcr_monarch] = RGB2YCbCr(rgb_monarch);

%% Displaying figures
figure('Name', 'Problem 2');
subplot(3,4,1); imshow(rgb_monarch); title('RGB');
subplot(3,4,2); imshow(R); title('Red');
subplot(3,4,3); imshow(G); title('Green');
subplot(3,4,4); imshow(B); title('Blue');

subplot(3,4,5); imshow(hsi_monarch); title('HSI');
subplot(3,4,6); imshow(H); title('Hue');
subplot(3,4,7); imshow(S); title('Saturation');
subplot(3,4,8); imshow(I); title('Intensity');

subplot(3,4,9); imshow(ycbcr_monarch); title('YCbCr');
subplot(3,4,10); imshow(Y); title('Y');
subplot(3,4,11); imshow(Cb); title('Cb');
subplot(3,4,12); imshow(Cr); title('Cr');
end

%%-------------------------------------------------------------%%
% Change RGB format of image
function [rA, gA, bA, A] = changeRGB(X)
[m, n] = size(X);
p = n/3;
A = zeros(m,p,3);

A(:,:,1) = X(:,1:3:(n-2));  % Red
A(:,:,2) = X(:,2:3:(n-1));  % Green
A(:,:,3) = X(:,3:3:n);      % Blue
A = uint8(A);

rA = zeros(m,p,3); rA(:,:,1) = A(:,:,1); rA = uint8(rA);
gA = zeros(m,p,3); gA(:,:,2) = A(:,:,2); gA = uint8(gA);
bA = zeros(m,p,3); bA(:,:,3) = A(:,:,3); bA = uint8(bA);
end

% Converting RGB to HSI
function [hA, sA, iA, A] = RGB2HSI(X)
[m, n, p] = size(X);
A = zeros(m, n, p);

X_norm = double(X)/255; % Normalize input RGB image

rX_norm = X_norm(:,:,1);    % Red normalized
gX_norm = X_norm(:,:,2);    % Green normalized
bX_norm = X_norm(:,:,3);    % Blue normalized

iA = (rX_norm + gX_norm + bX_norm)/3;   % Intensity
sA = 1 - min(X_norm,[],3)./iA;          % Saturation

% Calculate hue
V1 = 2 * rX_norm - gX_norm - bX_norm;
V2 = 2 * sqrt((rX_norm - gX_norm).^2 + (rX_norm - bX_norm).*(gX_norm - bX_norm));
theta = acos(V1./V2);
hA = theta/(2*pi);

A(:,:,1) = hA;
A(:,:,2) = sA;
A(:,:,3) = iA;
end

% Converting RGB to YCbCr
function [yA, cbA, crA, A] = RGB2YCbCr(X)
[m, n, p] = size(X);
A = zeros(m, n, p);

X_norm = double(X)/255; % Normalize input RGB image

for i = 1:m
   for j = 1:n
       % Y component
       A(i,j,1) = 0.299*X_norm(i,j,1) + 0.587*X_norm(i,j,2) + 0.114*X_norm(i,j,3);
       
       % Cb component
       A(i,j,2) = -0.169*X_norm(i,j,1) - 0.331*X_norm(i,j,2) + 0.500*X_norm(i,j,3) + 0.5;
       
       % Cr component
       A(i,j,3) = 0.500*X_norm(i,j,1) - 0.419*X_norm(i,j,2) - 0.081*X_norm(i,j,3) + 0.5;
   end
end

yA = A(:,:,1);
cbA = A(:,:,2);
crA = A(:,:,3);
end
