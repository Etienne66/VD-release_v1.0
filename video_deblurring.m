%% initialize
clear;
close all;
clc;
addpath('TV_L1_OF\');

img_dirs = ['dataset\car\input\'];
flow_dirs = ['dataset\car\flow\'];

n_imlist = 50;%number of imgs
imlist = dir([img_dirs '*.jpg']);%img extension
start_idx = 0;%start from first frame
pyramid_factor = 0.9;


%car%
lambda = 250;
tau = 0.23*ones(n_imlist,1);
pyramid_levels = 17;

% bicycle%
% lambda = 500;%1000;
% tau = 0.23*ones(n_imlist,1);
% pyramid_levels = 17;

%street%
% lambda = 500;
% tau = 0.23*ones(n_imlist,1);
% pyramid_levels = 17;

%bridge%
% lambda = 100%1000;
% tau = 0.32*ones(n_imlist,1);
% pyramid_levels = 17;

%books%
% lambda = 250;
% tau = 0.4*ones(n_imlist,1);
% pyramid_levels = 17;

%sign%
% lambda = 100;%%1000;%varying lambda from 25
% tau = 0.27*ones(n_imlist,1);
% pyramid_levels = 12;

%playground
% lambda = 250;
% tau = 0.35*ones(n_imlist,1);
% pyramid_levels = 14;



B = cell(n_imlist, 1);
d = cell(n_imlist, 1);
flow_forward = cell(n_imlist, 1);
flow_backward = cell(n_imlist, 1);


for i = 1 : n_imlist
    img_name = [img_dirs imlist(start_idx+i).name];
    B{i} = double(rgb2gray(imread(img_name)))/255.;
    B{i} = max(0, min(1, B{i}));
    
    if i < n_imlist
        flow_forward_name = [flow_dirs sprintf('%04dTo%04d.flo', start_idx+i, start_idx+i+1)];
        flow_forward{i} = readFlowFile(flow_forward_name);
    end
    if i > 1
        flow_backward_name = [flow_dirs sprintf('%04dTo%04d.flo',start_idx+i, start_idx+i-1)];
        flow_backward{i} = readFlowFile(flow_backward_name);
    end
end

flow_backward{1}(:,:,1) = -flow_forward{1}(:,:,1);
flow_backward{1}(:,:,2) = -flow_forward{1}(:,:,2);
flow_forward{n_imlist}(:,:,1) = -flow_backward{n_imlist}(:,:,1);
flow_forward{n_imlist}(:,:,2) = -flow_backward{n_imlist}(:,:,2);


% dual_temporal = zeros(M, N, n_imlist);
%% update
tic


err = cell(n_imlist, 1);
K = cell(n_imlist, 1);
L = B;
L_ = L;
% for i = 1 : n_imlist
%     if i < n_imlist
%         flow_forward_name = [flow_dirs sprintf('%04dTo%04d.flo', i, i+1)];
%         flow_forward{i} = readFlowFile(flow_forward_name);
%     end
%     if i > 1
%         flow_backward_name = [flow_dirs sprintf('%04dTo%04d.flo', i, i-1)];
%         flow_backward{i} = readFlowFile(flow_backward_name);
%     end
% end



[M, N, C] = size(L{1});
dual_spatial = zeros(M, N, n_imlist, 2);

num_temporal_neighbors = 2;
dual_temporal_forward = zeros(M, N, n_imlist, num_temporal_neighbors);
dual_temporal_backward = dual_temporal_forward;




%%




width_Pyramid = cell(pyramid_levels,1);
height_Pyramid = cell(pyramid_levels,1);

B_Pyramid = cell(n_imlist,1);


% precompute image sizes
width_Pyramid{1} = N;
height_Pyramid{1} = M;


for i = 2:pyramid_levels
    width_Pyramid{i} = pyramid_factor*width_Pyramid{i-1};
    height_Pyramid{i} = pyramid_factor*height_Pyramid{i-1};
    if min(width_Pyramid{i}, height_Pyramid{i}) < 16
        pyramid_levels = i;
        break;
    end
end

for i = 1:pyramid_levels
    width_Pyramid{i} = round(width_Pyramid{i});
    height_Pyramid{i} = round(height_Pyramid{i});
end


% set up image pyramides
for lv = pyramid_levels :-1:1;
    for n = 1 : n_imlist
        B_Pyramid{n} = imresize(B{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'bicubic');
        
        if lv == pyramid_levels
            L{n} = B_Pyramid{n};
            L_{n} = L{n};
            rescale_factor = width_Pyramid{pyramid_levels}/width_Pyramid{1};
            flow_forward{n} = rescale_factor*imresize(flow_forward{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'bicubic');
            flow_backward{n} = rescale_factor*imresize(flow_backward{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'bicubic');
        else
            if n == 1 || n == n_imlist
                L{n} = B_Pyramid{n};
            else
                L{n} = imresize(L{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'bicubic');
            end
            L_{n} = L{n};
            rescale_factor = width_Pyramid{lv}/width_Pyramid{lv+1};
            flow_forward{n} = rescale_factor*imresize(flow_forward{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'nearest');
            flow_backward{n} = rescale_factor*imresize(flow_backward{n}, [height_Pyramid{lv} width_Pyramid{lv}], 'nearest');
            
        end
    end
    
    [M N C] = size(L{1});
    
    dual_spatial = zeros(M, N, n_imlist, 2);
    dual_temporal_forward = zeros(M, N, n_imlist, num_temporal_neighbors);
    dual_temporal_backward = dual_temporal_forward;
    
    for kk = 1 : 1
        for n = 1 : n_imlist
            in_name = sprintf('in_%04d.png', n);
            imwrite(uint8(255*min(1, max(0, B_Pyramid{n}))), in_name);
        end
        
        for n = 2 : n_imlist - 1
            fprintf('%d ', n);
            K{n} = makeKernelMatrixFromFlows(flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau(n));            
            err{n} = detect_small_structure(L{n}, B_Pyramid{n}, K{n});            
        end
        
        [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = ...
            updateLatentImages(K, L, L_, B_Pyramid, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda, err);
        %         for iter = 1  : 5
        %             [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = ...
        %                 updateDuals2(K, L, L_, B_Pyramid, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda);
        %
        %             L = updatePrimals(L, K, B_Pyramid, num_temporal_neighbors, lambda, err);
        %             L_ = extrapolate(L, L_);
        %
        %             for n = 1 : n_imlist
        %                 out_name = sprintf('out1_%04d.png', n);
        %                 imwrite(uint8(255*min(1, max(0, L{n}))), out_name);
        %             end
        %
        %             fprintf('iter:%d\n', iter);
        %         end
        

            L = nonLocalMeans(L, flow_forward, flow_backward, 5, 3);
            for n = 1 : n_imlist
                out_name = sprintf('out2_%04d.png', n);                
                imwrite(uint8(255*min(1, max(0, L{n}))), out_name);                
            end

        
        
        
        [flow_forward, flow_backward] = updateFlows3(L, B_Pyramid, flow_forward, flow_backward, num_temporal_neighbors, tau);        
        
    end
    fprintf('Lv:%d\n', lv);
%     lambda = lambda/pyramid_factor
    
end

%% final color deblur
img_r = cell(n_imlist, 1);
img_g = cell(n_imlist, 1);
img_b = cell(n_imlist, 1);
for n = 1 : n_imlist
    img_name = [img_dirs imlist(start_idx+n).name];
    img = double(imread(img_name))./255.;
    img_r{n} = img(:,:,1);
    img_g{n} = img(:,:,2);
    img_b{n} = img(:,:,3);
end
for n = 2 : n_imlist - 1
    fprintf('%d ', n);
    K{n} = makeKernelMatrixFromFlows(flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau(n));
    err{n} = detect_small_structure(L{n}, B_Pyramid{n}, K{n});
end

%%

lambda = max(lambda, 250);
% L_r = cell(n_imlist,1);
% L_g = cell(n_imlist,1);
% L_b = cell(n_imlist,1);
% for n = 2 : n_imlist - 1
%     n
% L_r{n} = updateLatentImage(img_r{n}, img_r{n}, K{n}, lambda, 25);
% L_g{n} = updateLatentImage(img_g{n}, img_g{n}, K{n}, lambda, 25);
% L_b{n} = updateLatentImage(img_b{n}, img_b{n}, K{n}, lambda, 25);
% end
[L_r] = updateLatentImages(K, L, L_, img_r, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda, err);
[L_g] = updateLatentImages(K, L, L_, img_g, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda, err);
[L_b] = updateLatentImages(K, L, L_, img_b, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda, err);

% %%
L_r = nonLocalMeans(L_r, flow_forward, flow_backward, 3, 3);
L_g = nonLocalMeans(L_g, flow_forward, flow_backward, 3, 3);
L_b = nonLocalMeans(L_b, flow_forward, flow_backward, 3, 3);



for n = 1 : n_imlist
    img = zeros(M, N, 3);
    img(:,:,1) = L_r{n};
    img(:,:,2) = L_g{n};
    img(:,:,3) = L_b{n};
    out_name = sprintf('L_%04d.png', n);
    imwrite(uint8(255*min(1, max(0, img))), out_name);
end

%%


% TEST = nonLocalMeans(L, flow_forward, flow_backward, 3, 3);
% for n = 1 : n_imlist
%     out_name = sprintf('out2_%04d.png', n);
%     imwrite(uint8(255*min(1, max(0, TEST{n}))), out_name);
% end




