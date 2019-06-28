function [K]= makeKernelMatrixFromFlows(u_forward, v_forward,u_backward, v_backward, tau)

[M, N] = size(u_forward);
NUM_FRAME = 70;
[i, j, k] = mx_make_kernel_from_flows(u_forward, v_forward, u_backward, v_backward, tau, NUM_FRAME);
% M*N
% M*N*(1+2*NUM_FRAME)
% nnz(i)
% nnz(j)

% pause

% tic
K = sparse(i, j, 1./(1+2*NUM_FRAME), M*N, M*N);
% K = sparse(i, j, k, M*N, M*N);
% tic
% sum_K = sum(K');
% size(sum_K)
% sum_K = repmat(sum_K, N*M, 1);
% K = K./sum_K;
% toc


% tocT
% fprintf('1\n');
% 
% tic
% KL = mx_calc_KL(i,j,1./(1+2*NUM_FRAME), u_forward);
% toc
% fprintf('2\n');
% 
% tic
% K*reshape(u_forward', [M*N 1]);
% toc
% fprintf('3\n');

% K = sparse(i, j, 1, M*N, M*N);
% K = speye(M*N, M*N);

end


