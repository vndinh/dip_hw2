clear; clc;

Problem1();

function Problem1()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw2\Test_images');
file = fopen(fullfile(imgdir,'\Gray_monarch_noisy_512x512.raw'),'rb');
gray_monarch = fread(file,fliplr([512,512]),'*uint8')';
fclose(file);

%%
%%---------------------Insert code below ----------------------%%
H1 = [1 2 1; 2 4 2; 1 2 1]/16;
H2 = [1 4 1; 4 -20 4; 1 4 1]/16;
gray_monarch_1 = conv2Dsame(gray_monarch,H1);
gray_monarch_2 = conv2Dsame(gray_monarch,H2);

%% Displaying figures
figure('Name', 'Problem 1');
subplot(1,3,1); imshow(gray_monarch,[]); title('Gray Monarch original');
subplot(1,3,2); imshow(gray_monarch_1,[]); title('Kernel H1');
subplot(1,3,3); imshow(gray_monarch_2,[]); title('Kernel H2');
hold on;
end

%%-------------------------------------------------------------%%
% Convolution function that size of output image equals to input image
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