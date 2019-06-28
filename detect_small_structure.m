function [err, u, v] = detect_small_structure(L, B, K, delta)
if nargin ~= 4
    %     delta = 0.075;
    delta = 0.05;
%     delta = 0.035;
end



[M N C] = size(L);


rho = K*reshape(L', [M*N 1]);
rho = reshape(rho, [N M])' - B;
rho = delta - abs(rho);


L1 = 8;
tau =  1/sqrt(L1);
sigma = 1/tau/L1;

w = zeros(M, N);
w_ = w;

dual_w = zeros(M, N);
dual_w_x = zeros(M, N);
dual_w_y = zeros(M, N);
lambda = 50;% 50;%50;%250;
for iter = 1 : 50
    
    if 0
        dual_w = dual_w + sigma*w_;
        one_over_reprojection = 1./max(1.0, sqrt(dual_w.^2));
        dual_w = dual_w.*one_over_reprojection;
        
        w_ = w;
        div = dual_w;
    else
        w_x = dxp(w_);
        w_y = dyp(w_);
        
        dual_w_x = dual_w_x + sigma*w_x;
        dual_w_y =dual_w_y + sigma*w_y;
        
        one_over_reprojection = 1./max(1.0, sqrt(dual_w_x.^2 + dual_w_y.^2));
        dual_w_x = dual_w_x.*one_over_reprojection;
        dual_w_y = dual_w_y.*one_over_reprojection;
        
        w_ = w;
        div = dxm(dual_w_x) + dym(dual_w_y);
    end
    w = w + tau*(div - rho*lambda);
    
    
    w = max(0, min(1, w));
    w_ = 2*w - w_;
end


% err = zeros(M, N);
err = ones(M, N);
for i = 1 : M
    for j = 1 : N
        if w(i, j) > 0.5
            err(i, j) = 0.2;
%             err(i, j) = 0.01;
        end
    end
end





end


