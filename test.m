clear;
close all;


addpath(genpath('TV_L1_OF'));






img_dirs = ['C:\Papers\Video deblurring\video_deblur_dataset\bicycle\input\'];
flow_dirs = ['C:\Papers\Video deblurring\video_deblur_dataset\bicycle\flow\'];


imlist = dir([img_dirs '*.jpg']);
n_imlist = 4;%length(imlist)

n_imlist = 4
L = cell(n_imlist, 1);
d = cell(n_imlist, 1);
flow_forward = cell(n_imlist, 1);
flow_backward = cell(n_imlist, 1);

for i = 1 : n_imlist
    img_name = [img_dirs imlist(i).name];    
    L{i} = double(rgb2gray(imread(img_name)))/255.;
    
    if i < n_imlist
        flow_forward_name = [flow_dirs sprintf('%04dTo%04d.flo', i-1, i-1+1)];
        flow_forward{i} = readFlowFile(flow_forward_name);        
    end
    if i > 1
        flow_backward_name = [flow_dirs sprintf('%04dTo%04d.flo', i-1, i-1-1)];
        flow_backward{i} = readFlowFile(flow_backward_name);        
    end
end

L_diff_backward = mx_calc_temporal_diff(L{4}, L{3}, flow_backward{4}(:,:,1), flow_backward{4}(:,:,2));
figure
imshow(abs(L_diff_backward), [])
figure
imagesc(abs(L_diff_backward))
                    drawnow;
                    fprintf('kk\n');
%                     pause

    
% %     motion('0150.jpg', '0151.jpg', '150To151.flo');
%     
%     motion('0003.jpg', '0002.jpg', '3To2.flo');
%     
% 
% 
