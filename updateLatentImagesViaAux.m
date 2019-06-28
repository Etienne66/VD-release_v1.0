function [L] = updateLatentImagessViaAux(B, L, K, lambda, max_iter)
[M N] = size(B);
MN = M*N;

SIGMA = 10;%/10;%10%0.3536;%
TAU = 1/sqrt(8)^2/SIGMA;

[B_dx, B_dy] = extract_gradient(B);
% max_B_dx = max(abs(B_dx(:)));
% max_B_dy = max(abs(B_dy(:)));
% thr = 0.05;
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
e = zeros(M, N);
e_dx = zeros(M, N);
e_dy = zeros(M, N);
% e_= e;

dual_p_x = zeros(M, N);
dual_p_y = zeros(M, N);
% dual_e = zeros(M, N);
% dual_e_dx = zeros(M, N);
% dual_e_dy = zeros(M, N);

% if nargin < 6
%     g_x = ones(M, N);
%     g_y = ones(M, N);
% end
% w_I = 0.05;
w_I = 0;%0.25;

theta = 1./lambda;%0.01
one_over_2_theta = 1./(2*theta);
    two_tau_over_lambda = 2*TAU/lambda;

e_ = e;
e_dx_ = e_dx;
e_dy_ = e_dy;

for iter_out = 1 : 3

    
%     iter_out
    L_ = L;
    dual_e = SIGMA*e_;
    reprojection = max(1, sqrt(dual_e.^2));
    dual_e = dual_e./reprojection;
    
    dual_e_dx = SIGMA*e_dx_;
    reprojection = max(1, sqrt(dual_e_dx.^2));
    dual_e_dx = dual_e_dx./reprojection;
    
    dual_e_dy = SIGMA*e_dy_;
    reprojection = max(1, sqrt(dual_e_dy.^2));
    dual_e_dy = dual_e_dy./reprojection;
    
    
    [L_dx L_dy] = extract_gradient(L);    
    err = K*reshape(L', [MN 1]) - reshape(B', [MN 1]);
    err = reshape(err, [N M])';
    err_dx = K*reshape(L_dx', [MN 1]) - reshape(B_dx', [MN 1]);
    err_dx = reshape(err_dx, [N M])';
    err_dy = K*reshape(L_dy', [MN 1]) - reshape(B_dy', [MN 1]);
    err_dy = reshape(err_dy, [N M])';
    
    e_ = e;
    e_dx_ = e_dx;
    e_dy_ = e_dy;
    div_e = dual_e;
    div_e_dx = dual_e_dx;
    div_e_dy = dual_e_dy;
    
       
    e = e + TAU * div_e;
    e_dx = e_dx + TAU * div_e_dx;
    e_dy = e_dy + TAU * div_e_dy;
    
    e = (e - two_tau_over_lambda*err)/(1 -two_tau_over_lambda);
    e_dx = (e_dx - two_tau_over_lambda*err_dx)/(1 -two_tau_over_lambda);
    e_dy = (e_dy - two_tau_over_lambda*err_dy)/(1 -two_tau_over_lambda);
    
    e_ = 2*e - e_;
    e_dx_ = 2*e_dx - e_dx_;
    e_dy_ = 2*e_dy - e_dy_;
    
    
%     theta = theta/4.;
%     figure;
%     imagesc(e)

    for iter_in = 1 : max_iter
         
        %     iter
        u_x = dxp(L_);
        u_y = dyp(L_);
        
        dual_p_x = (dual_p_x + SIGMA*u_x);%/(1+SIGMA*epsilon);
        dual_p_y = (dual_p_y + SIGMA*u_y);%/(1+SIGMA*epsilon);
        
        
        if 1
            reprojection = max(1, sqrt(dual_p_x.^2 + dual_p_y.^2));
            dual_p_x = dual_p_x./reprojection;
            dual_p_y = dual_p_y./reprojection;
        else
            dual_p_x = max(-g_x, min(g_x, dual_p_x));
            dual_p_y = max(-g_y, min(g_y, dual_p_y));
        end
        
        div_u = dxm(dual_p_x) + dym(dual_p_y);
        L_ = L;
        L = L + TAU * div_u;
        [L_dx, L_dy] = extract_gradient(L);
        
        b = L/(2*TAU*one_over_2_theta) + (b_dx - e_dx) + (b_dy - e_dy) + w_I*(b_I - e);
        
        tmp_dx = K* reshape(L_dx', [N*M 1]);
        tmp_dx = reshape(K_t*tmp_dx, [N M])';
                
        tmp_dy = K* reshape(L_dy', [N*M 1]);
        tmp_dy = reshape(K_t*tmp_dy, [N M])';        
        
        tmp = K* reshape(L', [N*M 1]);
        tmp = reshape(K_t*tmp, [N M])';        
        
        AL = nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy ) + w_I*tmp;
        
        
        AL = (1/(2*TAU*one_over_2_theta))*L + AL;
        
        
        r = (b - AL);
        p = r;
        r_sqr = sum(sum(r.*r));
       
        tic;
        for(iter_inner = 1 : 5)
            
            [p_dx p_dy] = extract_gradient(p);
            
            tmp_dx = K* reshape(p_dx', [N*M 1]);
            tmp_dx = reshape(K_t*tmp_dx, [N M])';
            
            tmp_dy = K* reshape(p_dy', [N*M 1]);
            tmp_dy = reshape(K_t*tmp_dy, [N M])';
            
            tmp = K* reshape(p', [N*M 1]);
            tmp = reshape(K_t*tmp, [N M])';
            
            Ap = nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy ) + w_I*tmp;
            Ap = (1/(2*TAU*one_over_2_theta))*p + Ap;
            
            alpha = r_sqr/sum(sum(p.*Ap));
            
            r = r - alpha*Ap;
            
            tmp = sum(sum(r.*r));
            
            beta = 1/(r_sqr);
            r_sqr = tmp;
            
            L = L + alpha*p;
            
            beta = beta*r_sqr;
            p = r + beta*p;
        end
        toc
%         pause;
        L(L > 1) = 1;
    L(L < 0) = 0;
    
    L_ = 2*L - L_;
    
    L_(L_ > 1) = 1;
    L_(L_ < 0) = 0;
    
        
     
    end
    
       
        imshow(L_);
    tmp_title = sprintf('iter=%d', iter_in);
    title(tmp_title);
    drawnow;
end
%     L_ = bfilter2(L_, 5, [2.5 7/255]);
%     L = bfilter2(L_, 9, [4.5 5/255]);
%     L_ = L;

% imshow(L_);
% drawnow;

L = imsharpen(L_);
L = min(1, max(0, L));

% imwrite(uint8(L*255), 'L_.png');

end



