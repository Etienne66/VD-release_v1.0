function [L] = nonLocalMeans(L, flow_forward, flow_backward, block_size, window_size)

num_imgs = length(L);

% w = 2; % 5 x 5 window
% sigma = [2 10./255];
% sigma = [2 15./255];
% sigma = [3 20./255];

% for n = 2 : num_imgs - 1
%     B = L{n};
%     B(B > 1) = 1;
%     B(B < 0) = 0;
%     L{n} = bfilter2(B, w, sigma);
% 
% end

[M N C] = size(L{1});
for n = 2 : num_imgs - 1
% n = 17
    ref_img = L{n};
    qry_img_forward = L{n+1};
    qry_img_backward = L{n-1};
    
    [w0 v0] = mx_calc_non_local_means(ref_img, ref_img, zeros(M, N), zeros(M, N), block_size, window_size);
    [w_forward v_forward] = mx_calc_non_local_means(ref_img, qry_img_forward, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), block_size, window_size);
    [w_backward v_backward] = mx_calc_non_local_means(ref_img, qry_img_backward, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), block_size, window_size);

%     occ_forward = detectOcclusion_Brown(flow_forward{n}(:,:,1),flow_forward{n}(:,:,2));
%     occ_backward = detectOcclusion_Brown(flow_backward{n}(:,:,1),flow_backward{n}(:,:,2));
%     [w_forward, v_forward] = mx_calc_non_local_means_with_occ(ref_img, qry_img_forward, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), block_size, window_size, occ_forward);
%     [w_backward, v_backward] = mx_calc_non_local_means_with_occ(ref_img, qry_img_backward, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), block_size, window_size, occ_backward);
    
    if n == 2
        w = [w0; 0.25*w_forward];
        v = [v0; v_forward];
        sum_w = sum(w);
        sum_w = repmat(sum_w, 2*window_size*window_size, 1);
    elseif n == num_imgs - 1
        w = [w0; 0.25*w_backward];
        v = [v0; v_backward];
        sum_w = sum(w);
        sum_w = repmat(sum_w, 2*window_size*window_size, 1);
    else
    w = [w0; 0.25*w_forward; 0.25*w_backward];
    v = [v0; v_forward; v_backward];
    sum_w = sum(w);
    sum_w = repmat(sum_w, 3*window_size*window_size, 1);
    end
    w = w./sum_w;
    L{n} = reshape(sum(w.*v), [M N]);
    
    L{n} = imsharpen(L{n});
    L{n} = min(1, max(0, L{n}));
end






end



