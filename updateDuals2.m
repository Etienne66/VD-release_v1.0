function [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = updateDuals2(K, L, L_, B, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda )
mu = 2;
SIGMA = 10*num_temporal_neighbors;% 10;%200;%200;
TAU = 1/sqrt(8)^2/SIGMA/max(1,mu);

[M N C] = size(L{1});
num_imgs = length(L);

rho = cell(num_imgs, 1);
rho{1} = zeros(M, N);
rho{num_imgs} = zeros(M, N);
for n = 2 : num_imgs - 1
    tmp_rho = K{n}*reshape(L{n}', [M*N 1]);
    tmp_rho = abs(reshape(tmp_rho, [N M])' - B{n}).^0.8;    
    rho{n} = exp(-10*tmp_rho);
%     rho{n} = 1;
end


% for n = 1 : num_imgs
for n = 2 : num_imgs - 1
    L_x = dxp(L_{n});
    L_y = dyp(L_{n});

    dual_spatial(:,:,n,1) = (dual_spatial(:,:,n,1) + SIGMA*L_x);
    dual_spatial(:,:,n,2) = (dual_spatial(:,:,n,2) + SIGMA*L_y);
    
    reprojection = max(1, sqrt(dual_spatial(:,:,n,1).^2 + dual_spatial(:,:,n,2).^2));
    dual_spatial(:,:,n,1) = dual_spatial(:,:,n,1)./reprojection;
    dual_spatial(:,:,n,2) = dual_spatial(:,:,n,2)./reprojection;

    
    
      
    
    % compute divergence    
    spatial_div_L =  dxm(dual_spatial(:,:,n,1)) + dym(dual_spatial(:,:,n,2));    
    temporal_div_L = 0;
    
    u_forward = zeros(M, N);
    v_forward = zeros(M, N);
    u_backward = zeros(M, N);
    v_backward = zeros(M, N);
    
    m_forward = ones(M, N);
    m_backward = ones(M, N);
    
    for t = 1 : num_temporal_neighbors
        temporal_idx = n + t - 1;
        if temporal_idx < num_imgs
            [u_forward, v_forward, m_forward] = mx_calc_temporal_flow(u_forward, v_forward, flow_forward{temporal_idx}(:,:,1), flow_forward{temporal_idx}(:,:,2), m_forward);
            
            dual_temporal_forward(:,:,n,t) = dual_temporal_forward(:,:,n,t) + mu*SIGMA*mx_calc_temporal_diff(L_{n}, L_{n+t}, u_forward, v_forward);
%             reprojection = max(1, abs(dual_temporal_forward(:,:,n,t)));
%             dual_temporal_forward(:,:,n,t) = dual_temporal_forward(:,:,n,t)./reprojection;
            dual_temporal_forward(:,:,n,t) = max(-rho{n+t}, min(rho{n+t}, dual_temporal_forward(:,:,n,t)));
            
            idxx = repmat([1:N], M,1) + u_forward;
            idyy = repmat([1:M]',1,N) + v_forward;
            tmp_m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);
            tmp = dual_temporal_forward(:,:,n,t);
            tmp(tmp_m) = 0.0;
            tmp(~m_forward) = 0.0;
            dual_temporal_forward(:,:,n,t) = tmp;
            
            temporal_div_L = temporal_div_L + mx_calc_temporal_transposed_diff(dual_temporal_forward(:,:,n,t), u_forward, v_forward);     
        end
        
         temporal_idx = n - t + 1;
         if temporal_idx > 1
             [u_backward, v_backward, m_backward] = mx_calc_temporal_flow(u_backward, v_backward, flow_backward{temporal_idx}(:,:,1), flow_backward{temporal_idx}(:,:,2), m_backward);
             
             dual_temporal_backward(:,:,n,t) = dual_temporal_backward(:,:,n,t) + mu*SIGMA*mx_calc_temporal_diff(L_{n}, L_{n-t}, u_backward, v_backward);
             dual_temporal_backward(:,:,n,t) = max(-rho{n-t}, min(rho{n-t}, dual_temporal_backward(:,:,n,t)));
             
             idxx = repmat([1:N], M,1) + u_backward;
             idyy = repmat([1:M]',1,N) + v_backward;
             tmp_m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);
             tmp = dual_temporal_backward(:,:,n,t);
             tmp(tmp_m) = 0.0;
             tmp(~m_backward) = 0.0;
             dual_temporal_backward(:,:,n,t) = tmp;
             
             temporal_div_L = temporal_div_L + mx_calc_temporal_transposed_diff(dual_temporal_backward(:,:,n,t), u_backward, v_backward);
         end
    end    
    
    % remember old L
    L_{n} = L{n};
    L{n} = L{n} + TAU*(spatial_div_L - mu*temporal_div_L);
end

end



