%%
clear;
close all;

addpath(genpath('C:\Papers\Video deblurring\video_deblur_dataset\street\input'));
addpath(genpath('C:\Papers\Video deblurring\video_deblur_dataset\street\flow'));
addpath(genpath('TV_L1_OF'));
addpath(genpath('fastdeconv1'));

img = double(imread('0005.jpg'))./255;
[M N C] = size(img);




flow_150_149 = readFlowFile('0005to0004.flo');
flow_150_151 = readFlowFile('0005to0006.flo');
[M N C] = size(flow_150_151)


u_forward = flow_150_151(:,:,1);
v_forward = flow_150_151(:,:,2);
u_backward = flow_150_149(:,:,1);
v_backward = flow_150_149(:,:,2);
% tau = 0.23;% 0.25;%0.23;
tau = 0.23;% 0.25;%0.23;



tic
K = makeKernelMatrixFromFlows(u_forward, v_forward, u_backward, v_backward, tau);
toc



fprintf('done');
%%


fast_hyper_Laplacian_deconv(K, img)


% 
% tic
% out = img;
% [M N C] = size(img);
% tic
% for c = 1 : C
%         out(:,:,c) = updateLatentImageViaAux(img(:,:,c), img(:,:,c), K, 1000, 5);
% %     out(:,:,c) = updateLatentImage(img(:,:,c), img(:,:,c), K, 100, 5);
% %     out(:,:,c) = updateLatentImage2(img(:,:,c), img(:,:,c), u_forward, v_forward, u_backward, v_backward, tau, 100, 5);
% end
% toc
% fprintf('done\n');
% imshow(uint8(255*out));
% imwrite(uint8(255*out), 'out.png');








