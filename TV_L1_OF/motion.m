function [] = motion(img1_name, img2_name, flo_name)

% img_src1 = getappdata(handles.figure_motion,'img_src1');
% img_src2 = getappdata(handles.figure_motion,'img_src2');
img_src1 = imread(img1_name);
img_src2 = imread(img2_name);

[M N C] = size(img_src1);
if C == 3
    img_src1 = double(rgb2gray(img_src1))/255.;
    img_src2 = double(rgb2gray(img_src2))/255.;
end

lambda = 40;
beta   = 0.01;
max_iter = 50

pyramid_levels = 1000;
pyramid_factor = 0.9;
warps = 1;

[flow illumination] = coarse_to_fine(img_src1, img_src2, lambda, beta, warps, max_iter, pyramid_levels, pyramid_factor);



tmp = flow;
% find robust max flow for better visualization
magnitude = (tmp(:,:,1).^2 + tmp(:,:,2).^2).^0.5;  
max_flow = prctile(magnitude(:),95);


tmp(:,:,1) = min(max(tmp(:,:,1),-max_flow),max_flow);
tmp(:,:,2) = min(max(tmp(:,:,2),-max_flow),max_flow);

imwrite(uint8(flowToColor(tmp)),'motion.png');
% illumination = illumination - min(illumination(:));
% illumination = illumination/max(illumination(:));
% imwrite(illumination,'illumination.png');
writeFlowFile(tmp, flo_name);
% append the estimated optical flow to the first frame to show the motion
% of every pixel



end
