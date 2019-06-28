clear;
close all;

addpath(genpath('TV_L1_OF'));

img_dirs = ['C:\Papers\Video deblurring\video_deblur_dataset\bridge\input\'];
flow_dirs = ['C:\Papers\Video deblurring\video_deblur_dataset\bridge\flow\'];

imlist = dir([img_dirs '*.jpg']);
n_imlist = length(imlist)



% flow_forward = cell(n_imlist, 1);
% flow_backward = cell(n_imlist, 1);

% for i = 1 : n_imlist
i = 82;
flow_forward_name = [flow_dirs sprintf('%04dTo%04d.flo', i, i+1)];
flow_backward_name = [flow_dirs sprintf('%04dTo%04d.flo', i, i-1)];
flow_forward = readFlowFile(flow_forward_name);
flow_backward = readFlowFile(flow_backward_name);


[M N C] = size(flow_backward)

kernel = double(imread('kernel.png'));
% neg_dix = kernel < 0.05*max(kernel(:));
% kernel(neg_dix) = 0;
kernel = kernel./sum(kernel(:));
[kernel_size, ~] = size(kernel);


num_samples = 25;
one_over_num_samples = 1./num_samples;
one_over_1plus_2_times_num_samples =  1./(2*num_samples + 1);
kernel_default = zeros(kernel_size, kernel_size);
half_kernel_size = ceil(kernel_size/2);
kernel_default(half_kernel_size, half_kernel_size) = one_over_1plus_2_times_num_samples;

% tau = 0.5
for tau = 0.05 : 0.005: 0.5
    %     tau = 1
    
    kernel_sum = zeros(kernel_size, kernel_size);
    %     for i = 421 : 5 : 620 %bicycle
    %         for j = 161 : 5 : 360
    for i = 181 : 5 : 380 % car
        for j = 421 : 5 : 620
            
            
            kernel = kernel_default;
            
            for s = 1 : num_samples
                
                u = tau*flow_backward(i, j, 1)*s*one_over_num_samples;
                v = tau*flow_backward(i, j, 2)*s*one_over_num_samples;
                
                u = max(1, min(kernel_size, round(half_kernel_size + u)));
                v = max(1, min(kernel_size, round(half_kernel_size + v)));
                
                kernel(u, v) = kernel(u, v) + one_over_1plus_2_times_num_samples;
                
                u = tau*flow_forward(i, j, 1)*s*one_over_num_samples;
                v = tau*flow_forward(i, j, 2)*s*one_over_num_samples;
                
                u = max(1, min(kernel_size, round(half_kernel_size + u)));
                v = max(1, min(kernel_size, round(half_kernel_size + v)));
                
                kernel(u, v) = kernel(u, v) + one_over_1plus_2_times_num_samples;
            end
            %
            
            kernel_sum = kernel_sum + kernel;
        end
    end
    g = fspecial('gaussian', [3 3], 0.3);
    kernel_sum = imfilter(kernel_sum, g, 'same');
    
    
    
    
    max_kernel = max(kernel_sum(:));
    neg_idx = 0.05*max_kernel > max_kernel;
    kernel_sum(neg_idx) = 0;
    kernel_sum = kernel_sum./sum(kernel_sum(:));
    dist = abs(kernel_sum - kernel);
    dist = dist.*dist;
    
    ham_dist = 0;
    for kx = 1 : kernel_size
        for ky = 1 : kernel_size
            if kernel(ky, kx) > 0 && kernel_sum(ky, kx) > 0
            elseif kernel(ky, kx) == 0 && kernel_sum(ky, kx) == 0
            else
                ham_dist = ham_dist + 1;
            end
        end
    end
    
    
    kernel_sum = 255*kernel_sum./max(kernel(:));
    imwrite(kernel_sum, 'kernel_sum.png')
    
    drawnow;
    fprintf('tau:%f, dist:%f ham_dist:%f\n', tau, sum(dist(:)), ham_dist);
%     pause
end



%
%


