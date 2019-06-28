function [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = updateDuals(K, L, L_, B, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda )

SIGMA = 10*num_temporal_neighbors;% 10;%200;%200;
TAU = 1/sqrt(8)^2/SIGMA;

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
    
    if 1
    reprojection = max(1, sqrt(dual_spatial(:,:,n,1).^2 + dual_spatial(:,:,n,2).^2));
    dual_spatial(:,:,n,1) = dual_spatial(:,:,n,1)./reprojection;
    dual_spatial(:,:,n,2) = dual_spatial(:,:,n,2)./reprojection;
    else
        thr = lambda/1000.;
        dual_spatial(:,:,n,1) = max(-thr, min(thr, dual_spatial(:,:,n,1)));
        dual_spatial(:,:,n,2) = max(-thr, min(thr, dual_spatial(:,:,n,1)));
    end
    
    
      
    
    % compute divergence    
    spatial_div_L =  dxm(dual_spatial(:,:,n,1)) + dym(dual_spatial(:,:,n,2));    
    temporal_div_L = 0;
    
    
    for t = 1 : num_temporal_neighbors
        
        
        u_forward = zeros(M, N);
        v_forward = zeros(M, N);
        u_backward = zeros(M, N);
        v_backward = zeros(M, N);
        
        m_forward = ones(M, N);
        m_backward = ones(M, N);        
     
        
        for tt = 1 : t
            temporal_idx = n + tt - 1;
            if temporal_idx < num_imgs
%             if temporal_idx < num_imgs - 1
                [u_forward, v_forward, m_forward] = mx_calc_temporal_flow(u_forward, v_forward, flow_forward{temporal_idx}(:,:,1), flow_forward{temporal_idx}(:,:,2), m_forward);
            else
                m_forward = [];
            end
%             fprintf('%f %f\n', u_forward(261, 4),v_forward(261, 4));
%             pause;
            temporal_idx = n - tt + 1;
            if temporal_idx > 1
%             if temporal_idx > 2
                [u_backward, v_backward, m_backward] = mx_calc_temporal_flow(u_backward, v_backward, flow_backward{temporal_idx}(:,:,1), flow_backward{temporal_idx}(:,:,2), m_backward);
            else
                m_backward = [];
            end
            
        end
        %
        
        
        if~isempty(m_forward)
            
            idxx = repmat([1:N], M,1) + u_forward;
            idyy = repmat([1:M]',1,N) + v_forward;
            tmp_m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);
            tmp_m = ~tmp_m;
            
            for i = 1 : M
                for j = 1 : N
                    if tmp_m(i,j) == 0 || m_forward(i, j) == 0
                        m_forward(i, j) = 0;
                    end
                end
            end
            
            dual_temporal_forward(:,:,n,t) = dual_temporal_forward(:,:,n,t) + SIGMA*mx_calc_temporal_diff(L_{n}, L_{n+t}, u_forward, v_forward);
%             reprojection = max(1, abs(dual_temporal_forward(:,:,n,t)));
%             dual_temporal_forward(:,:,n,t) = dual_temporal_forward(:,:,n,t)./reprojection;
            dual_temporal_forward(:,:,n,t) = max(-rho{n+t}, min(rho{n+t}, dual_temporal_forward(:,:,n,t)));
  
            
            dual_temporal_forward(:,:,n,t) = mx_handle_boundary(dual_temporal_forward(:,:,n,t), m_forward);
            temporal_div_L = temporal_div_L + mx_calc_temporal_transposed_diff(dual_temporal_forward(:,:,n,t), u_forward, v_forward);            
        end
        if~isempty(m_backward)
            idxx = repmat([1:N], M,1) + u_backward;
            idyy = repmat([1:M]',1,N) + v_backward;
            tmp_m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);
            tmp_m = ~tmp_m;
            
            for i = 1 : M
                for j = 1 : N
                    if tmp_m(i,j) == 0 || m_backward(i, j) == 0
                        m_backward(i, j) = 0;
                    end
                end
            end
            
            dual_temporal_backward(:,:,n,t) = dual_temporal_backward(:,:,n,t) + SIGMA*mx_calc_temporal_diff(L_{n}, L_{n-t}, u_backward, v_backward);
%             reprojection = max(1, abs(dual_temporal_backward(:,:,n,t)));
%             dual_temporal_backward(:,:,n,t) = dual_temporal_backward(:,:,n,t)./reprojection;
            dual_temporal_backward(:,:,n,t) = max(-rho{n-t}, min(rho{n-t}, dual_temporal_backward(:,:,n,t)));
            
            dual_temporal_backward(:,:,n,t) = mx_handle_boundary(dual_temporal_backward(:,:,n,t), m_backward);
            temporal_div_L = temporal_div_L + mx_calc_temporal_transposed_diff(dual_temporal_backward(:,:,n,t), u_backward, v_backward);
        end
    end
    
    % remember old L
    L_{n} = L{n};
    L{n} = L{n} + TAU*(spatial_div_L - temporal_div_L);
end

end



