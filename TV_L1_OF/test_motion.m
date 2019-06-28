clear;
close all;


addpath('..\deblur');

img_dirs = ['..\dataset\car\input\'];
flow_dirs = ['..\dataset\car\flow\'];
imlist = dir([img_dirs '*.png']);
n_imlist = length(imlist)

img_name = [img_dirs imlist(1).name];
img_forward_name = [img_dirs imlist(2).name];    
motion(img_name, img_forward_name, [flow_dirs sprintf('%04dTo%04d.flo', 1, 2)]);
    
for i = 2 : n_imlist-1    
    img_name = [img_dirs imlist(i).name];    
    
    if i < n_imlist
        img_forward_name = [img_dirs imlist(i+1).name];
        motion(img_name, img_forward_name, [flow_dirs sprintf('%04dTo%04d.flo', i, i+1)]);
    end
    
    if i > 1
        img_backward_name = [img_dirs imlist(i-1).name];
        motion(img_name, img_backward_name, [flow_dirs sprintf('%04dTo%04d.flo', i, i-1)]);
    end
    
end

img_name = [img_dirs imlist(n_imlist).name];
img_backward_name = [img_dirs imlist(n_imlist-1).name];    
motion(img_name, img_backward_name, [flow_dirs sprintf('%04dTo%04d.flo', n_imlist, n_imlist-1)]);