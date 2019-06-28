%%
clear;
close all;

addpath(genpath('TV_L1_OF'));


flow_150_149 = readFlowFile('flow_150_149.flo');
flow_150_151 = readFlowFile('flow_150_151.flo');
[M N C] = size(flow_150_149)

u_150_151 = flow_150_151(:,:,1);
v_150_151 = flow_150_151(:,:,2);
u_150_149 = flow_150_149(:,:,1);
v_150_149 = flow_150_149(:,:,2);

% 
% 0.23*u_150_151(582,340)
% 0.23*v_150_151(582,340)
% 0.23*u_150_149(582,340)
% 0.23*v_150_149(582,340)
% pause
% 
% u_150_151 = 5*ones(M,N);
% v_150_151 = 3*ones(M,N);
% u_150_149 = -3*ones(M,N);
% v_150_149 = -2*ones(M,N);







tau = 0.23;
% tau = 0.05;
% [M, N] = size(u);
% histogram = zeros(M*N, 1);

numFrame = 100;%70;
num_coeff = (1+2*numFrame);

x_idx = zeros(M*N*num_coeff, 1);
y_idx = zeros(M*N*num_coeff, 1);
z_idx = zeros(M*N*num_coeff, 1);
one_over_num_frame = 1./(1+2*numFrame);
count = 1;

for j = 1:M
    for i = 1:N
        k = (j-1)*N + i;
        jj = j;
        ii = i;
        kk = (jj - 1)*N + ii;
        x_idx(count) = k;
        y_idx(count) = kk;
        z_idx(count) = one_over_num_frame;
        count = count + 1;
        
        for f = 1:numFrame
            
            jj = max(1, min(M, round(j + tau*v_150_151(j, i)*f/numFrame)));
            ii = max(1, min(N, round(i + tau*u_150_151(j, i)*f/numFrame)));
            kk = (jj - 1)*N + ii;
            x_idx(count) = k;
            y_idx(count) = kk;
            z_idx(count) = one_over_num_frame;
            count = count + 1;
            
            

            jj = max(1, min(M, round(j + tau*v_150_149(j, i)*f/numFrame)));
            ii = max(1, min(N, round(i + tau*u_150_149(j, i)*f/numFrame)));
            kk = (jj - 1)*N + ii;
            x_idx(count) = k;
            y_idx(count) = kk;
            z_idx(count) = one_over_num_frame;
            count = count + 1;
        end
        
          
    end
end
count = count - 1;

K = sparse(x_idx(1:count), y_idx(1:count), 1./num_coeff, M*N, M*N);
%%

img = double(imread('0150.jpg'))./255;
out = img;
[M N C] = size(img);
for c = 1 : C
%     out(:,:,c) = updateLatentImageViaAux(img(:,:,c), img(:,:,c), K, 500, 50)
out(:,:,c) = updateLatentImage(img(:,:,c), img(:,:,c), K, 1000, 30);
end

fprintf('done');
%%
imshow(uint8(255*out));
imwrite(uint8(255*out), 'out.png');






