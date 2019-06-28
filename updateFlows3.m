function [flow_forward, flow_backward] = updateFlows3(L, B, flow_forward, flow_backward, num_temporal_neighbors, exposure_time)
%
% tau = 1/sqrt(8);
% sigma = 1/sqrt(8);

sigma = 10;
tau = 1/sqrt(8)^2/sigma;

% sigma = 1;
% tau = 1/sqrt(8)^2/sigma;

num_imgs = length(L);

% lambda = 50/4.;%3*2.5;
lambda = 50/4.;% means nu = 4*lambda / 50
% lambda = 7.5;%1;%7.5;%7.5;% means nu = 4*lambda / 50
% lambda = 50;
max_iter = 20;%15;


fixed_flow_forward = flow_forward;
fixed_flow_backward = flow_backward;

[M N C] = size(L{1});

for n = 1 : num_imgs
    fprintf('OF:%d\n', n);
    [L_dx, L_dy] = extract_gradient(L{n});
        g_x = ones(M, N);
        g_y = ones(M, N);
        var = 2*(25/255)^2;
        for i = 1 : M - 1
            for j = 1 : N - 1
                g_x(i,j) = max(0.001, exp(-L_dx(i,j)*L_dx(i,j)/var));
                g_y(i,j) = max(0.001, exp(-L_dy(i,j)*L_dy(i,j)/var));
            end
        end
    
    %     B1 = B{n};
    [B_dx, B_dy] = extract_gradient(B{n});
    
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
                u_accum = u;
                v_accum = v;
                m_accum = ones(M, N);
                
                for t = 1 : num_temporal_neighbors
                    temporal_idx = n + t - 1;
                    if temporal_idx < num_imgs
                        if t > 1
                            [u_accum, v_accum, m_accum] = mx_calc_temporal_flow(u_accum, v_accum, fixed_flow_forward{temporal_idx}(:,:,1), fixed_flow_forward{temporal_idx}(:,:,2), m_accum);
                        end
                        
                        L2 = L{n+t};
                        rho_du_p = rho_du_p + calc_rho(L{n}, L2, u_accum+du, v_accum, m_accum);
                        rho_du_m = rho_du_m + calc_rho(L{n}, L2, u_accum-du, v_accum, m_accum);
                        rho_dv_p = rho_dv_p + calc_rho(L{n}, L2, u_accum, v_accum+dv, m_accum);
                        rho_dv_m = rho_dv_m + calc_rho(L{n}, L2, u_accum, v_accum-dv, m_accum);
                    end
                end
                
                %
                %                [KL] = mx_calc_KL(L{n}, u + du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
                %                rho_du_p = rho_du_p + abs(KL - B{n});
                %
                %                [KL] = mx_calc_KL(L{n}, u - du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
                %                rho_du_m = rho_du_m + abs(KL - B{n});
                %
                %                [KL] = mx_calc_KL(L{n}, u, v + dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
                %                rho_dv_p = rho_du_p + abs(KL - B{n});
                %
                %                [KL] = mx_calc_KL(L{n}, u, v - dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), tau);
                %                rho_dv_m = rho_dv_m + abs(KL - B{n});
                
                ww = 1;
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, u + du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), exposure_time(n));
                rho_du_p = rho_du_p + ww*(0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2);
                
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, u - du, v, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), exposure_time(n));
                rho_du_m = rho_du_m + ww*(0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2);
                %
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, u, v + dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), exposure_time(n));
                rho_dv_p = rho_dv_p + ww*(0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2);
                
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, u, v - dv, flow_backward{n}(:,:,1), flow_backward{n}(:,:,2), exposure_time(n));
                rho_dv_m = rho_dv_m + ww*(0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2);
                
            end
            if k == 2
                u_accum = u;
                v_accum = v;
                m_accum = ones(M, N);
                
                for t = 1 : num_temporal_neighbors
                    temporal_idx = n - t + 1;
                    if temporal_idx > 1
                        if t > 1
                            [u_accum, v_accum, m_accum] = mx_calc_temporal_flow(u_accum, v_accum, fixed_flow_backward{temporal_idx}(:,:,1), fixed_flow_backward{temporal_idx}(:,:,2), m_accum);
                        end
                        
                        L2 = L{n-t};
                        rho_du_p = rho_du_p + calc_rho(L{n}, L2, u_accum+du, v_accum, m_accum);
                        rho_du_m = rho_du_m + calc_rho(L{n}, L2, u_accum-du, v_accum, m_accum);
                        rho_dv_p = rho_dv_p + calc_rho(L{n}, L2, u_accum, v_accum+dv, m_accum);
                        rho_dv_m = rho_dv_m + calc_rho(L{n}, L2, u_accum, v_accum-dv, m_accum);
                    end
                end
                
                %
                %                 [KL] = mx_calc_KL(L{n}, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u + du, v, tau);
                %                 rho_du_p = rho_du_p + abs(KL - B{n});
                %
                %                 [KL] = mx_calc_KL(L{n}, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u - du, v, tau);
                %                 rho_du_m = rho_du_m + abs(KL - B{n});
                %
                %                 [KL] = mx_calc_KL(L{n}, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v + dv, tau);
                %                 rho_dv_p = rho_dv_p + abs(KL - B{n});
                %
                %                 [KL] = mx_calc_KL(L{n}, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v - dv, tau);
                %                 rho_dv_m = rho_dv_m + abs(KL - B{n});
                
                %
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u + du, v, exposure_time);
                rho_du_p = rho_du_p + 0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2;
                
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u - du, v, exposure_time);
                rho_du_m = rho_du_m + 0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2;
                
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v + dv, exposure_time);
                rho_dv_p = rho_dv_p + 0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2;
                
                [KL, KL_dx, KL_dy] = mx_calc_KL_KLdx_KLdy(L{n}, L_dx, L_dy, flow_forward{n}(:,:,1), flow_forward{n}(:,:,2), u, v - dv, exposure_time);
                rho_dv_m = rho_dv_m + 0.25*abs(KL - B{n}).^2 + abs(KL_dx - B_dx).^2 + abs(KL_dy - B_dy).^2;
                
                
                
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
                p(:,:,1) = (p(:,:,1) + sigma*g_x.*u_x);
                p(:,:,2) = (p(:,:,2) + sigma*g_y.*u_y);
                
                p(:,:,3) = (p(:,:,3) + sigma*g_x.*v_x);
                p(:,:,4) = (p(:,:,4) + sigma*g_y.*v_y);
                
                %             % reprojection to |pu| <= 1
                reprojection_uv = max(1.0, sqrt(p(:,:,1).^2 + p(:,:,2).^2 + p(:,:,3).^2 + p(:,:,4).^2));
                p(:,:,1) = p(:,:,1)./reprojection_uv;
                p(:,:,2) = p(:,:,2)./reprojection_uv;
                p(:,:,3) = p(:,:,3)./reprojection_uv;
                p(:,:,4) = p(:,:,4)./reprojection_uv;
                %
                %             p(:,:,1) = max(-g_x, min(g_x, p(:,:,1)));
                %             p(:,:,2) = max(-g_y, min(g_y, p(:,:,2)));
                %             p(:,:,3) = max(-g_x, min(g_x, p(:,:,3)));
                %             p(:,:,4) = max(-g_y, min(g_y, p(:,:,4)));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PRIMAL
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % remember old u,v,w
                u_ = u;
                v_ = v;
                %             w_ = w;
                
                % compute divergence
                div_u = dxm(g_x.*p(:,:,1)) + dym(g_y.*p(:,:,2));
                div_v = dxm(g_x.*p(:,:,3)) + dym(g_y.*p(:,:,4));
                %             div_w = dxm(p(:,:,5)) + dym(p(:,:,6));
                
                % update u,v,w
                u = u + tau*(div_u);
                v = v + tau*(div_v);
                
                u = u - tau*lambda.*rho_du;
                v = v - tau*lambda.*rho_dv;
                
                u_ = 2*u - u_;
                v_ = 2*v - v_;
                
                
                %             if mod(iter,15) == 0
%                 if iter == 15
%                     show_flow(u,v,L{n});
%                     fprintf('nth frame:%d tv-l1-motion-primal-dual: it = %d\n', n, iter)
%                 end
                
            end
%             u = medfilt2(u, [5 5], 'symmetric');
%             v = medfilt2(v, [5 5], 'symmetric');            
        end
        
      
%         occ = detectOcclusion_Brown(u, v);
%         uv = cat(length(size(u))+1, u, v);
%         uv = denoise_color_weighted_medfilt2(uv, 255*L{n}, occ, 7, 5, 15, 0);
%         u = uv(:,:,1);
%         v = uv(:,:,2);
        
        u = medfilt2(u, [5 5], 'symmetric');
        v = medfilt2(v, [5 5], 'symmetric');

        
        if k == 1 && n < num_imgs
            flow_forward{n}(:,:,1) = u;
            flow_forward{n}(:,:,2) = v;
        elseif k == 2 && n > 1
            flow_backward{n}(:,:,1) = u;
            flow_backward{n}(:,:,2) = v;
        end
        
    end
    
    if n < num_imgs
        imwrite(uint8(flowToColor(flow_forward{n})),sprintf('M%04dto%04d.png', n, n+1));
        writeFlowFile(flow_forward{n}, sprintf('M%04dto%04d.flo', n, n+1));
    end
    if n > 1
        imwrite(uint8(flowToColor(flow_backward{n})),sprintf('M%04dTo%04d.png', n, n-1));
        writeFlowFile(flow_backward{n}, sprintf('M%04dTo%04d.flo', n, n-1));
    end
    


end

end



