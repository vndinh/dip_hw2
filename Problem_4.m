clear; clc;

Problem4();

function Problem4()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw2\Test_images');
file = fopen(fullfile(imgdir,'\Color_lenna_256x256.raw'),'rb');
color_lenna = fread(file,fliplr([256,256*3]),'*uint8')';
fclose(file);

%%
%%---------------------Insert code below ----------------------%%
% changing data format to RGB image
[~, ~, ~, lenna_256x256_original] = changeRGB(color_lenna);

% Kernel for 4:1 magnification
Peg = [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]; % Peg kernel
Pyramid = conv2D(Peg,Peg)/16;               % Pyramid kernel
Cubic = conv2D(Pyramid,Pyramid)/16;         % Cubic B-spline kernel

% LPF before down sampling
lenna_256x256_lpfD_peg = lpf(lenna_256x256_original,Peg/16);
lenna_256x256_lpfD_pyramid = lpf(lenna_256x256_original,Pyramid/16);
lenna_256x256_lpfD_cubic = lpf(lenna_256x256_original,Cubic/16);

% Down sampling
lenna_down64x64_peg = downSampling(lenna_256x256_lpfD_peg,Peg,4,4);
lenna_down64x64_pyramid = downSampling(lenna_256x256_lpfD_pyramid,Pyramid,4,4);
lenna_down64x64_cubic = downSampling(lenna_256x256_lpfD_cubic,Cubic,4,4);

% Up sampling
lenna_up256x256_peg = upSampling(lenna_down64x64_peg,4,4);
lenna_up256x256_pyramid = upSampling(lenna_down64x64_pyramid,4,4);
lenna_up256x256_cubic = upSampling(lenna_down64x64_cubic,4,4);

% LPF after up sampling
lenna_256x256_lpfU_peg = lpf(lenna_up256x256_peg,Peg);
lenna_256x256_lpfU_pyramid = lpf(lenna_up256x256_pyramid,Pyramid);
lenna_256x256_lpfU_cubic = lpf(lenna_up256x256_cubic,Cubic);

psnr_peg = PSNR(lenna_256x256_original,lenna_256x256_lpfU_peg);
psnr_pyramid = PSNR(lenna_256x256_original,lenna_256x256_lpfU_pyramid);
psnr_cubic = PSNR(lenna_256x256_original,lenna_256x256_lpfU_cubic);

%% Displaying result
figure('Name', 'Problem 4');
subplot(3,3,2); imshow(uint8(lenna_256x256_original)); title('Original Lenna');
subplot(3,3,4); imshow(uint8(lenna_down64x64_peg)); title('Down Sampling Peg');
subplot(3,3,5); imshow(uint8(lenna_down64x64_pyramid)); title('Down Sampling Pyramid');
subplot(3,3,6); imshow(uint8(lenna_down64x64_cubic)); title('Down Sampling Cubic B-spline');
subplot(3,3,7); imshow(uint8(lenna_256x256_lpfU_peg)); title('Up Sampling Peg');
subplot(3,3,8); imshow(uint8(lenna_256x256_lpfU_pyramid)); title('Up Sampling Pyramid');
subplot(3,3,9); imshow(uint8(lenna_256x256_lpfU_cubic)); title('Up Sampling Cubic B-spline');

fprintf('PSNR_Peg = %f dB\n',psnr_peg);
fprintf('PSNR_Pyramid = %f dB\n',psnr_pyramid);
fprintf('PSNR_Cubic = %f dB\n',psnr_cubic);

end

%%-------------------------------------------------------------%%
% Separating RGB components
function [rA, gA, bA, A] = changeRGB(X)
[m, n] = size(X);
p = n/3;

rA = zeros(m,p,3); rA(:,:,1) = X(:,1:3:(n-2)); rA = uint8(rA);  % Red
gA = zeros(m,p,3); gA(:,:,2) = X(:,2:3:(n-1)); gA = uint8(gA);  % Green
bA = zeros(m,p,3); bA(:,:,3) = X(:,3:3:n); bA = uint8(bA);      % Blue

A = rA + gA + bA;    % RGB image
end

% Low pass filter
function A = lpf(X,H)
[xRows, xCols, xColor] = size(X);
A = zeros(xRows,xCols,xColor);

% RGB seperating
xR = zeros(xRows,xCols); xR(:,:) = X(:,:,1);
xG = zeros(xRows,xCols); xG(:,:) = X(:,:,2);
xB = zeros(xRows,xCols); xB(:,:) = X(:,:,3);

% Low pass filter
A(:,:,1) = conv2Dsame(xR,H);
A(:,:,2) = conv2Dsame(xG,H);
A(:,:,3) = conv2Dsame(xB,H);
end

function A = downSampling(X,H,m,n)
[xRows, xCols, xColor] = size(X);

aRows = floor(xRows/m);
aCols = floor(xCols/n);
A = zeros(aRows, aCols, xColor);

A(:,:,:) = X(1:m:xRows, 1:n:xCols, :);
end

function A = upSampling(X,m,n)
[xRows, xCols, xColor] = size(X);

aRows = xRows * m;
aCols = xCols * n;
A = zeros(aRows,aCols,xColor);

A(1:m:aRows, 1:n:aCols, :) = X(:,:,:);
end

% Convolution function that size of output image equals to of input image
function A = conv2Dsame(X,H)
[xRows, xCols] = size(X);
[hRows, hCols] = size(H);

% Rotate kernel matrix by 180 degree
K = rot90(H,2);

% Zero padding matrix
center = floor((size(K)+1)/2);
left = center(2) - 1;
right = hCols - center(2);
top = center(1) - 1;
bot = hRows - center(1);
Z = zeros(xRows + top + bot, xCols + left + right);

for i = (1+top) : (xRows+top)
    for j = (1+left) : (xCols+left)
        Z(i,j) = X(i - top, j - left);
    end
end

% Convolution
A = zeros(xRows,xCols);
for m = 1:xRows
   for n = 1:xCols
      for i = 1:hRows
         for j = 1:hCols
            p = m - 1;
            q = n - 1;
            A(m,n) = A(m,n) + Z(i + p,j + q) * K(i,j);
         end
      end
   end
end
end

% Convolution function that size of output image is greater than of input image
function Y = conv2D(X,H)
[xR, xC] = size(X);
[hR, hC] = size(H);

% Rotate kernel matrix by 180 degree
K = rot90(H,2);

% Size of ouput matrix
yR = xR + hR - 1;
yC = xC + hC - 1;
Y = zeros(yR,yC);

% Zero padding matrix
zR = yR + hR;
zC = yC + hC;
Z = zeros(zR,zC);
for i = hR:(hR + xR - 1)
    for j = hC:(hC + xC - 1)
        Z(i,j) = X(i - hR + 1, j - hC + 1);
    end
end

% Convolution
for m = 1:yR
   for n = 1:yC
      for i = 1:hR
         for j = 1:hC
            p = m - 1;
            q = n - 1;
            Y(m,n) = Y(m,n) + Z(i + p,j + q) * K(i,j);
         end
      end
   end
end
end

function psnr = PSNR(Orig,Dist)
[m, n, p] = size(Orig);
Orig = double(Orig);
Dist = double(Dist);
error = Orig - Dist;
MSE = sum(sum(sum(error.^2)))/(m*n*p);
if MSE > 0
    psnr = 20*log10(max(max(max(Orig)))) - 10*log10(MSE);
else
    psnr = 99;
end
end

