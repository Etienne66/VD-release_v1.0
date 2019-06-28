function [flow_backward] = updateBackwardFlows(L, flow_backward)

tau = 1/sqrt(8);
sigma = 1/sqrt(8);

num_imgs = length(L);

lambda = 40;
beta   = 0.01;
max_iter = 50;



for n = 2 : num_imgs
    [M N C] = size(L{n});    
    
    I1 = L{n};%double(rgb2gray(img_src1))/255.;
    I2 = L{n-1};%double(rgb2gray(img_src2))/255.;
    
    u = flow_backward{n}(:,:,1);
    v = flow_backward{n}(:,:,2);    
    
    
%     g_x = ones(M, N);
%     g_y = ones(M, N);
%     exp_var = exp(-1./((7/255.)^2));
%     for i = 1 : M-1
%         for j = 1 : N-1
%             dist_x = (I1(i, j) - I1(i, j+1))^2;
%             g_x(i, j) = exp(dist_x)*exp_var;
%             
%             dist_y = (I1(i, j) - I1(i+1, j))^2;
%             g_y(i, j) = exp(dist_y)*exp_var;
%         end
%     end
    
    
    % vectorization
    u_ = u;
    v_ = v;
    w = zeros(M, N);
    w_ = w;
    p = zeros(M, N, 6);
     
    u0 = u;
    v0 = v;
    
    [I_x, I_y, I_t, I_2_warped] = warping(I1, I2, u0, v0);
    I_grad_sqr = max(1e-09, I_x.^2 + I_y.^2 + beta*beta);
    
    
    for k = 0:max_iter-1        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DUAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute derivatives and update dual variable
        u_x = dxp(u_);
        u_y = dyp(u_);
        
        v_x = dxp(v_);
        v_y = dyp(v_);
        
        w_x = dxp(w_);
        w_y = dyp(w_);
        
        % update dual variable
        p(:,:,1) = (p(:,:,1) + sigma*u_x);
        p(:,:,2) = (p(:,:,2) + sigma*u_y);
        
        p(:,:,3) = (p(:,:,3) + sigma*v_x);
        p(:,:,4) = (p(:,:,4) + sigma*v_y);


        
        p(:,:,5) = (p(:,:,5) + sigma*w_x);
        p(:,:,6) = (p(:,:,6) + sigma*w_y);
        
        % reprojection to |pu| <= 1
        reprojection_uv = max(1.0, sqrt(p(:,:,1).^2 + p(:,:,2).^2 + p(:,:,3).^2 + p(:,:,4).^2));
        p(:,:,1) = p(:,:,1)./reprojection_uv;
        p(:,:,2) = p(:,:,2)./reprojection_uv;
        p(:,:,3) = p(:,:,3)./reprojection_uv;
        p(:,:,4) = p(:,:,4)./reprojection_uv;

%  p(:,:,1) = max(-g_x, min(g_x, p(:,:,1)));
%         p(:,:,2) = max(-g_y, min(g_y, p(:,:,2)));
%         p(:,:,3) = max(-g_x, min(g_x, p(:,:,3)));
%         p(:,:,4) = max(-g_y, min(g_y, p(:,:,4)));
        
        
        reprojection_w = max(1.0, sqrt(p(:,:,5).^2 + p(:,:,6).^2));
        p(:,:,5) = p(:,:,5)./reprojection_w;
        p(:,:,6) = p(:,:,6)./reprojection_w;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PRIMAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remember old u,v,w
        u_ = u;
        v_ = v;
        w_ = w;
        
        % compute divergence
        div_u = dxm(p(:,:,1)) + dym(p(:,:,2));
        div_v = dxm(p(:,:,3)) + dym(p(:,:,4));
        div_w = dxm(p(:,:,5)) + dym(p(:,:,6));
        
        % update u,v,w
        u = u + tau*(div_u);
        v = v + tau*(div_v);
        w = w + tau*(div_w);
        
        % prox operator for u,v,w
        rho = I_t + (u - u0).*I_x + (v - v0).*I_y + beta*w;
        idx1 = rho      < - tau*lambda*I_grad_sqr;
        idx2 = rho      >   tau*lambda*I_grad_sqr;
        idx3 = abs(rho) <=  tau*lambda*I_grad_sqr;
        
        u(idx1) = u(idx1) + tau*lambda*I_x(idx1);
        v(idx1) = v(idx1) + tau*lambda*I_y(idx1);
        w(idx1) = w(idx1) + tau*lambda*beta;
        
        u(idx2) = u(idx2) - tau*lambda*I_x(idx2);
        v(idx2) = v(idx2) - tau*lambda*I_y(idx2);
        w(idx2) = w(idx2) - tau*lambda*beta;
        
        u(idx3) = u(idx3) - rho(idx3).*I_x(idx3)./I_grad_sqr(idx3);
        v(idx3) = v(idx3) - rho(idx3).*I_y(idx3)./I_grad_sqr(idx3);
        w(idx3) = w(idx3) - rho(idx3).*beta./I_grad_sqr(idx3);
        
        u_ = 2*u - u_;
        v_ = 2*v - v_;
        w_ = 2*w - w_;
        
        if mod(k,10) == 0
            show_flow(u,v,I1,I2,I_2_warped + (u-u0).*I_x + (v-v0).*I_y + beta*w);
            fprintf('nth frame:%d tv-l1-motion-primal-dual: it = %d\n', n, k)
        end
        
    end
    
    
    
    
    %         occ = detectOcclusion_Brown(u, v);
    %
    %         uv = cat(length(size(u))+1, u, v);
    %         uv = denoise_color_weighted_medfilt2(uv, 255*L{n}, occ, 7, 5, 10, 0);
    %         u = uv(:,:,1);
    %         v = uv(:,:,2);
    
    flow_backward{n}(:,:,1) = medfilt2(u, [5 5], 'symmetric');
    flow_backward{n}(:,:,2) = medfilt2(v, [5 5], 'symmetric');
end
end



