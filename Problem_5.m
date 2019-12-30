clear; clc;

Problem5();

function Problem5()
% Template for EE535 Digial Image Processing
% Insert the code in the designated area below
%% Loading directory for image files
imgdir = uigetdir('D:\KAIST\Courses\dip\hw\hw2\Test_images');
file = fopen(fullfile(imgdir,'\Gray_monarch_512x512.raw'),'rb');
gray_monarch = fread(file,fliplr([512,512]),'*uint8')';
fclose(file);

%%
%%---------------------Insert code below ----------------------%%
rand_noise = round(32*rand(512,512)-16*ones(512,512));
gray_monarch = double(gray_monarch);

% Uniform quantize 2 bit
un_quan_2 = uniform_quantizer(gray_monarch+rand_noise,2);
un_quan_4 = uniform_quantizer(gray_monarch+rand_noise,4);
un_quan_6 = uniform_quantizer(gray_monarch+rand_noise,6);

un_quan_2 = un_quan_2 - rand_noise;
un_quan_4 = un_quan_4 - rand_noise;
un_quan_6 = un_quan_6 - rand_noise;

% PSNR for uniform quantizier
un_psnr_2 = PSNR(gray_monarch,un_quan_2);
un_psnr_4 = PSNR(gray_monarch,un_quan_4);
un_psnr_6 = PSNR(gray_monarch,un_quan_6);

%% Displaying figures
figure('Name', 'Problem 5');
subplot(3,3,2); imshow(uint8(gray_monarch)); title('Original Gray Monarch');
subplot(3,3,4); imshow(uint8(un_quan_2)); title('Uniform Quantizer 2 bit');
subplot(3,3,5); imshow(uint8(un_quan_4)); title('Uniform Quantizer 4 bit');
subplot(3,3,6); imshow(uint8(un_quan_6)); title('Uniform Quantizer 6 bit');

fprintf('Uniform Quantizer:\n');
fprintf('2 bit: PSNR = %f dB\n', un_psnr_2);
fprintf('4 bit: PSNR = %f dB\n', un_psnr_4);
fprintf('6 bit: PSNR = %f dB\n', un_psnr_6);

%%---------------------------------------------------------------%%
end

% Compandor function
function A = uniform_quantizer(X, bit)  % input: 2D
    [m,n] = size(X);
    range = 255;
    L = 2^bit;
    q = range/L;
    t = 0:q:255;        % Decision(transition) levels
    r = (q/2):q:255;    % Reconsturction levels
    A = zeros(m,n);
    X = floor(X/q);     % Range: 0(t(1)),1,2,...,L-1(t(L)),L(t(L+1))
    
    for i = 0:length(t)-2
        [row, col]=find(X==i);
        for j = 1:length(row)
            A(row(j),col(j)) = r(i+1);
        end
    end
    
    [row, col] = find(X==length(t)-1);
    for j = 1:length(row)
            A(row(j),col(j)) = r(length(t)-1);
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
