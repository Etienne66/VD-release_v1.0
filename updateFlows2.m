function [flow_forward, flow_backward] = updateFlows2(L, B, flow_forward, flow_backward, tau, s_r)
% 
% tau = 1/sqrt(8);
% sigma = 1/sqrt(8);

sigma = 10;
tau = 1/sqrt(8)^2/sigma;

% sigma = 1;
% tau = 1/sqrt(8)^2/sigma;

num_imgs = length(L);

lambda = 50;
max_iter = 15;


[M N C] = size(L{1}); 

for n = 1 : num_imgs    
    L1 = L{n};
    [L1_dx, L1_dy] = extract_gradient(L1);
    
    B1 = B{n};
    [B1_dx, B1_dy] = extract_gradient(B1);
    
%     [P P_x P_y] = prediction(L1, 0.01,  0.1);
%     [L1, L1_dx, L1_dy] = prediction(L1, s_r,  0.1);
    
    for k = 1 : 2
        
        if k == 1 && n < num_imgs
            u = flow_forward{n}(:,:,1);
            v = flow_forward{n}(:,:,2);
        elseif k == 2 && n > 1            
            u = flow_backward{n}(:,:,1);
            v = flow_backward{n}(:,:,2);
        else
            continue;
        end
        
  
        % vectorization
        u_ = u;
        v_ = v;
        p = zeros(M, N, 4);
        
        du = 0.5;
        dv = 0.5;
            
        for warps = 1 : 5
            
            
        
            rho_du_p = 0;
            rho_du_m = 0;
            rho_dv_p = 0;
            rho_dv_m = 0;
            
            if k == 1
               
                
                if n < num_imgs
                    L2 = L{n+1};
                    [L2_dx, L2_dy] = extract_gradient(L2);                    
                   
                    rho_du_p = calc_rho(L1, L2, u+du, v);% + calc_rho(L1_dx, L2_dx, u+du, v) + calc_rho(L1_dy, L2_dy, u+du, v);
                    rho_du_m = calc_rho(L1, L2, u-du, v);% + calc_rho(L1_dx, L2_dx, u-du, v) + calc_rho(L1_dy, L2_dy, u-du, v);
                    rho_dv_p = calc_rho(L1, L2, u, v+dv);% + calc_rho(L1_dx, L2_dx, u, v+dv) + calc_rho(L1_dy, L2_dy, u, v+dv);
                    rho_dv_m = calc_rho(L1, L2, u, v-dv);% + calc_rho(L1_dx, L2_dx, u, v-dv) + calc_rho(L1_dy, L2_dy, u, v-dv);
                end
               
%                [KL] = mx_calc_KL(L1, u + du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                rho_du_p = rho_du_p + 0.8*abs(KL - B1);
%                
%                [KL] = mx_calc_KL(L1, u - du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                rho_du_m = rho_du_m + 0.8*abs(KL - B1);
%                
%                [KL] = mx_calc_KL(L1, u, v + dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                rho_dv_p = rho_du_p + 0.8*abs(KL - B1);
%                
%                [KL] = mx_calc_KL(L1, u, v - dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                rho_dv_m = rho_dv_m + 0.8*abs(KL - B1);
               
               
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, u + du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                 rho_du_p = rho_du_p + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
%                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, u - du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                 rho_du_m = rho_du_m + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
% %                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, u, v + dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                 rho_dv_p = rho_dv_p + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
%                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, u, v - dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
%                 rho_dv_m = rho_dv_m + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
                
            end
            if k == 2
                if n > 1
                    L2 = L{n-1};
                    [L2_dx, L2_dy] = extract_gradient(L2);                   

                    rho_du_p = calc_rho(L1, L2, u+du, v);% + calc_rho(L1_dx, L2_dx, u+du, v) + calc_rho(L1_dy, L2_dy, u+du, v);
                    rho_du_m = calc_rho(L1, L2, u-du, v);% + calc_rho(L1_dx, L2_dx, u-du, v) + calc_rho(L1_dy, L2_dy, u-du, v);
                    rho_dv_p = calc_rho(L1, L2, u, v+dv);% + calc_rho(L1_dx, L2_dx, u, v+dv) + calc_rho(L1_dy, L2_dy, u, v+dv);
                    rho_dv_m = calc_rho(L1, L2, u, v-dv);% + calc_rho(L1_dx, L2_dx, u, v-dv) + calc_rho(L1_dy, L2_dy, u, v-dv);
                end
                
%                
%                 [KL] = mx_calc_KL(L1, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u + du, v, tau);
%                 rho_du_p = rho_du_p + 0.8*abs(KL - B1);
%                 
%                 [KL] = mx_calc_KL(L1, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u - du, v, tau);
%                 rho_du_m = rho_du_m + 0.8*abs(KL - B1);
%                 
%                 [KL] = mx_calc_KL(L1, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v + dv, tau);
%                 rho_dv_p = rho_dv_p + 0.8*abs(KL - B1);
%                 
%                 [KL] = mx_calc_KL(L1, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v - dv, tau);
%                 rho_dv_m = rho_dv_m + 0.8*abs(KL - B1);
                
                
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u + du, v, tau);
%                 rho_du_p = rho_du_p + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
%                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u - du, v, tau);
%                 rho_du_m = rho_du_m + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
%                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v + dv, tau);
%                 rho_dv_p = rho_dv_p + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
%                 
%                 [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L1, L1_dx, L1_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v - dv, tau);
%                 rho_dv_m = rho_dv_m + abs(KL - B1) + abs(KL_dx - B1_dx) + abs(KL_dy - B1_dy);
                
                
                
            end
            
            rho_du = rho_du_p - rho_du_m;
            rho_dv = rho_dv_p - rho_dv_m;
        
        for iter = 0 : max_iter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DUAL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute derivatives and update dual variable
            u_x = dxp(u_);
            u_y = dyp(u_);
            
            v_x = dxp(v_);
            v_y = dyp(v_);
            
            % update dual variable
            p(:,:,1) = (p(:,:,1) + sigma*u_x);
            p(:,:,2) = (p(:,:,2) + sigma*u_y);
            
            p(:,:,3) = (p(:,:,3) + sigma*v_x);
            p(:,:,4) = (p(:,:,4) + sigma*v_y);
                        
            % reprojection to |pu| <= 1
            reprojection_uv = max(1.0, sqrt(p(:,:,1).^2 + p(:,:,2).^2 + p(:,:,3).^2 + p(:,:,4).^2));
            p(:,:,1) = p(:,:,1)./reprojection_uv;
            p(:,:,2) = p(:,:,2)./reprojection_uv;
            p(:,:,3) = p(:,:,3)./reprojection_uv;
            p(:,:,4) = p(:,:,4)./reprojection_uv;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PRIMAL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remember old u,v,w
            u_ = u;
            v_ = v;
%             w_ = w;
            
            % compute divergence
            div_u = dxm(p(:,:,1)) + dym(p(:,:,2));
            div_v = dxm(p(:,:,3)) + dym(p(:,:,4));
%             div_w = dxm(p(:,:,5)) + dym(p(:,:,6));
            
            % update u,v,w
            u = u + tau*(div_u);
            v = v + tau*(div_v);
            
            u = u - tau*lambda.*rho_du;
            v = v - tau*lambda.*rho_dv;
            
            u_ = 2*u - u_;
            v_ = 2*v - v_;
            
            
            if mod(iter,15) == 0
                show_flow(u,v,L{n});
                fprintf('nth frame:%d tv-l1-motion-primal-dual: it = %d\n', n, iter)
            end
            
        end
            u = medfilt2(u, [5 5], 'symmetric');
            v = medfilt2(v, [5 5], 'symmetric');
        end
        
        if k == 1 && n < num_imgs
            flow_forward{n}(:,:,1) = u;
            flow_forward{n}(:,:,2) = v;
        elseif k == 2 && n > 1
            flow_backward{n}(:,:,1) = u;
            flow_backward{n}(:,:,2) = v;
        end

    end
end

end



