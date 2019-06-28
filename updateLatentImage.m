function [L] = updateLatentImage(B, L, K, lambda, max_iter)
[M N] = size(B);

epsilon = 0;%7./255;

SIGMA = 10;%/10;%10%0.3536;%
TAU = 1/sqrt(8)^2/SIGMA;

[B_dx, B_dy] = extract_gradient(B);
% max_B_dx = max(abs(B_dx(:)));
% max_B_dy = max(abs(B_dy(:)));
% thr = 0.05;
% % thr = 0.1;
% B_dx(abs(B_dx) < thr*max_B_dx) = 0;
% B_dy(abs(B_dy) < thr*max_B_dy) = 0;

K_t = K';


tmp_Kt_B = K_t*reshape(B', [M*N 1]);
tmp_Kt_Bdx = K_t*reshape(B_dx', [M*N 1]);
tmp_Kt_Bdy = K_t*reshape(B_dy', [M*N 1]);

b_I = ( reshape(tmp_Kt_B, [N M])' );
b_dx = nablax_t_A( reshape(tmp_Kt_Bdx, [N M])' );
b_dy = nablay_t_A( reshape(tmp_Kt_Bdy, [N M])' );
% end

L_ = L;
dual_p_x = zeros(M, N);
dual_p_y = zeros(M, N);



% w_I = 0.05;%0.01;
w_I = 0.25;
for(iter = 1 : max_iter)
    %     iter
    u_x = dxp(L_);
    u_y = dyp(L_);
    
    dual_p_x = (dual_p_x + SIGMA*u_x);%/(1+SIGMA*epsilon);
    dual_p_y = (dual_p_y + SIGMA*u_y);%/(1+SIGMA*epsilon);
    
    if 1
        reprojection = max(1, sqrt(dual_p_x.^2 + dual_p_y.^2));
        dual_p_x = dual_p_x./reprojection;
        dual_p_y = dual_p_y./reprojection;  
    end
    
    div_u = dxm(dual_p_x) + dym(dual_p_y);
    L_ = L;
    L = L + TAU * div_u;
    [L_dx, L_dy] = extract_gradient(L);
    
    b = L/(2*TAU*lambda) + b_dx + b_dy + w_I*b_I;
        
    tmp_dx = K* reshape(L_dx', [N*M 1]);
    tmp_dx = reshape(K_t*tmp_dx, [N M])';
    
    tmp_dy = K* reshape(L_dy', [N*M 1]);
    tmp_dy = reshape(K_t*tmp_dy, [N M])';
    
    tmp = K* reshape(L', [N*M 1]);
    tmp = reshape(K_t*tmp, [N M])';
    
    AL = nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy ) + w_I*tmp;
    
    AL = (1/(2*TAU*lambda))*L + AL;
    
    
    r = (b - AL);
    p = r;
    r_sqr = sum(sum(r.*r));
    
    for(iter_inner = 1 : 5)
        
        [p_dx p_dy] = extract_gradient(p);
            
        tmp_dx = K* reshape(p_dx', [N*M 1]);        
        tmp_dx = reshape(K_t*tmp_dx, [N M])';
        
        tmp_dy = K* reshape(p_dy', [N*M 1]);        
        tmp_dy = reshape(K_t*tmp_dy, [N M])';
        
        tmp = K* reshape(p', [N*M 1]);        
        tmp = reshape(K_t*tmp, [N M])';
        
        Ap = nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy ) + w_I*tmp;
        Ap = (1/(2*TAU*lambda))*p + Ap;
        
        alpha = r_sqr/sum(sum(p.*Ap));
        
        r = r - alpha*Ap;
        
        tmp = sum(sum(r.*r));
        
        beta = 1/(r_sqr);
        r_sqr = tmp;
        
        L = L + alpha*p;
        
        beta = beta*r_sqr;
        p = r + beta*p;       
    end
    
    L(L > 1) = 1;
    L(L < 0) = 0;
    
    L_ = 2*L - L_;
    
%     L_(L_ > 1) = 1;
%     L_(L_ < 0) = 0;
    

end


 L = imsharpen(L);
 L = min(1, max(0, L));
        
    imshow(L_);
    tmp_title = sprintf('iter=%d', iter);
    title(tmp_title);
    drawnow;

end



