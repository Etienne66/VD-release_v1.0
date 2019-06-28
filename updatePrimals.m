function [set_L] = updatePrimals(set_L, set_K, set_B, num_temporal_neighbors, lambda, set_err)

SIGMA = 10*num_temporal_neighbors;%10*num_temporal_neighbors;;
TAU = 1/sqrt(8)^2/SIGMA;

[M N C] = size(set_L{1});
num_imgs = length(set_L);

w_I = 0.25;
% w_I = 0.025;
for n = 2 : num_imgs - 1
    K = set_K{n};
    K_t = K';
    L = set_L{n};
    B = set_B{n};
    err = set_err{n};
    err = reshape(err', [M*N 1]);
    
    [B_dx, B_dy] = extract_gradient(B);

    tmp_Kt_B = K_t*(err.*reshape(B', [M*N 1]));
    tmp_Kt_Bdx = K_t*(err.*reshape(B_dx', [M*N 1]));
    tmp_Kt_Bdy = K_t*(err.*reshape(B_dy', [M*N 1]));
    
   
    
%     tmp_Kt_B = tmp_Kt_B.*err;
%     tmp_Kt_Bdx = tmp_Kt_Bdx.*err;
%     tmp_Kt_Bdy = tmp_Kt_Bdy.*err;
    
    tmp_Kt_B = reshape(tmp_Kt_B, [N M])';
    tmp_Kt_Bdx = reshape(tmp_Kt_Bdx, [N M])';
    tmp_Kt_Bdy = reshape(tmp_Kt_Bdy, [N M])';

    b_I = ( tmp_Kt_B );
    b_dx = nablax_t_A( tmp_Kt_Bdx );
    b_dy = nablay_t_A( tmp_Kt_Bdy );


    [L_dx, L_dy] = extract_gradient(L);    
    b = L/(2*TAU*lambda) + w_I*b_I + b_dx + b_dy;
    
    tmp = K* reshape(L', [N*M 1]);
    tmp = tmp.*err;
    tmp = reshape(K_t*tmp, [N M])';
    
        
    tmp_dx = K* reshape(L_dx', [N*M 1]);
    tmp_dx = tmp_dx.*err;
    tmp_dx = reshape(K_t*tmp_dx, [N M])';
    
    tmp_dy = K* reshape(L_dy', [N*M 1]);
    tmp_dy = tmp_dy.*err;
    tmp_dy = reshape(K_t*tmp_dy, [N M])';
    
    
    
    
    
    
    AL =  w_I*tmp + nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy );
    
    AL = (1/(2*TAU*lambda))*L + AL;
    
    
    r = (b - AL);
    p = r;
    r_sqr = sum(sum(r.*r));
    
    for(iter_inner = 1 : 5)
        
        [p_dx p_dy] = extract_gradient(p);
            
        tmp_dx = K* reshape(p_dx', [N*M 1]);        
        tmp_dx = tmp_dx.*err;
        tmp_dx = reshape(K_t*tmp_dx, [N M])';
        
        tmp_dy = K* reshape(p_dy', [N*M 1]);        
        tmp_dy = tmp_dy.*err;
        tmp_dy = reshape(K_t*tmp_dy, [N M])';
        
        tmp = K* reshape(p', [N*M 1]);  
        tmp = tmp.*err;
        tmp = reshape(K_t*tmp, [N M])';
        
        
        
        
        
        Ap = w_I*tmp + nablax_t_A(tmp_dx) + nablay_t_A( tmp_dy );
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
    
    set_L{n} = L;
end
    

end


