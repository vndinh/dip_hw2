clear; clc;

Problem3();

% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
function Problem3()
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw2\Test_images');
file = fopen(fullfile(imgdir,'\Color_baboon_256x256.raw'),'rb');
color_baboon = fread(file,fliplr([256,256*3]),'*uint8')';
fclose(file);

%%
%%---------------------Insert code below ----------------------%%
[~,~,~,baboon_original] = changeRGB(color_baboon);

% Prameter Setting
% Assume that (1,1) is the origin point
% x axis be j and y axis be i.
Tr_x = 30;  
Tr_y = -40;
Sc_x = 1.2;     % x axis (j)
Sc_y = 0.8;     % y axis (i)
Ang = pi/6;     % radian

baboon_translation = translation(baboon_original,Tr_x,Tr_y);
baboon_rotation = rotation(baboon_original,Ang);
baboon_scaling = scaling(baboon_original,Sc_x,Sc_y);

% Translation -> Rotation -> Scaling
baboon_T = translation(baboon_original,Tr_x,Tr_y);
baboon_TR = rotation(baboon_T,Ang);
baboon_TRS = scaling(baboon_TR,Sc_x,Sc_y);

%% Displaying figures
figure('Name','Problem 3');
subplot(2,3,1); imshow(baboon_original); title('Original baboon');
subplot(2,3,2); imshow(uint8(baboon_translation)); title('Translating baboon');
subplot(2,3,3); imshow(uint8(baboon_rotation)); title('Rotating baboon');
subplot(2,3,4); imshow(uint8(baboon_scaling)); title('Scaling baboon');
subplot(2,3,5); imshow(uint8(baboon_TRS)); title('T -> R -> S');
end

%%---------------------------------------------------------------%%
% Change RGB format of image
function [rA, gA, bA, A] = changeRGB(X)
[m, n] = size(X);
p = n/3;

rA = zeros(m,p,3); rA(:,:,1) = X(:,1:3:(n-2)); rA = uint8(rA);  % Red
gA = zeros(m,p,3); gA(:,:,2) = X(:,2:3:(n-1)); gA = uint8(gA);  % Green
bA = zeros(m,p,3); bA(:,:,3) = X(:,3:3:n); bA = uint8(bA);      % Blue

A = rA + gA + bA;    % RGB image
end

% Translation
function A = translation(X,Tr_x,Tr_y)
T = [1 0 Tr_x; 0 1 Tr_y; 0 0 1];    % Translation matrix
[xRows, xCols, xColor] = size(X);

% Size of output image
aRows = ceil(xRows+abs(Tr_x));
aCols = ceil(xCols+abs(Tr_y));
A = zeros(aRows,aCols,xColor);

for a_i = 1:aRows
    for a_j = 1:aCols
        % Pre-translated index using inverse T
        tmp_idx = T\[a_i,a_j,1]';
        x_i = round(tmp_idx(1));
        x_j = round(tmp_idx(2));
        if x_i>=1 && x_i<=xRows && x_j>=1 && x_j<=xCols
           A(a_i,a_j,:) = X(x_i,x_j,:); 
        end
    end
end
end

% Rotation
function A = rotation(X,Ang)
% Rotation matrix
Cosin = cos(Ang);
Sin = sin(Ang);
R = [Cosin Sin  0; -Sin Cosin 0; 0 0 1];

[xRows, xCols, xColor] = size(X);

% Size of output image
V = [1 xRows 1 xRows; 1 1 xCols xCols; 1 1 1 1];    % Vertex
rotV = R*V;

min_vx = min(rotV(1,:)); max_vx = max(rotV(1,:));
aRows = ceil(max_vx - min_vx);

min_vy = min(rotV(2,:)); max_vy = max(rotV(2,:));
aCols = ceil(max_vy - min_vy);

A = zeros(aRows, aCols, xColor);

shift_x = round(V(1,3));
for a_i = 1:aRows
   for a_j = 1:aCols
      % Pre-rotated index using inverse R
      tmp_idx = R\[a_i,a_j-shift_x,1]';
      x_i = round(tmp_idx(1));
      x_j = round(tmp_idx(2));
      if x_i>=1 && x_i<=xRows && x_j>=1 && x_j<=xCols
         A(a_i,a_j,:) = X(x_i,x_j,:);
      end
   end
end
end

% Scaling
function A = scaling(X,Sc_x,Sc_y)
[xRows, xCols, xColor] = size(X);

% Size of output image
aRows = ceil(xRows*Sc_x);
aCols = ceil(xCols*Sc_y);
A = zeros(aRows,aCols,xColor);

% Scaling matrix
S = [Sc_x 0 0; 0 Sc_y 0; 0 0 1];

if Sc_x<=0 || Sc_y<=0
    print('Error: Scaling parameters must be positive');
end

for a_i = 1:aRows
    for a_j = 1:aCols
       tmp_idx = S\[a_i,a_j,1]'; 
       x_i = round(tmp_idx(1));
       x_j = round(tmp_idx(2));
       if x_i>=1 && x_i<=xRows && x_j>=1 && x_j<=xCols
           A(a_i,a_j,:) = X(x_i,x_j,:);
       end
    end
end
end
