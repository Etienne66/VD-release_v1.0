%%
clear;
close all;

addpath(genpath('TV_L1_OF'));

img = double(imread('0150.png'))./255;
[M N C] = size(img);

kernel = double(imread('kernel2.png'));

for ky = 1 : 19
    for kx = 1 : 19
        dist = (ky - 9)^2 + (kx - 9)^2;
        
        if kernel(ky, kx) > 0
            kernel(ky, kx) = 1;%exp(-dist./25)            
        end
    end
end
kernel = kernel/sum(kernel(:))
imshow(kernel, [])
drawnow;
% pause
K = uniformKernel2SparseMat(kernel, M, N);
fprintf('done');
%%


out = img;
[M N C] = size(img);
for c = 1 : C
    %     out(:,:,c) = updateLatentImageViaAux(img(:,:,c), img(:,:,c), K, 500, 50)
    out(:,:,c) = updateLatentImage(img(:,:,c), img(:,:,c), K, 1000, 30);
end

fprintf('done\n');
%%
imshow(uint8(255*out));
imwrite(uint8(255*out), 'out.png');






