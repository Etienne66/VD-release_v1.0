function [K histogram]= makeKernelMatrixFromFlow(u, v)
%
[M, N] = size(u);
histogram = zeros(M*N, 1);

numFrame = 70;
% numFrame = 50;
% num_coeff = (1+2*numFrame*2);
num_coeff = (1+2*numFrame);

x_idx = zeros(M*N*num_coeff, 1);
y_idx = zeros(M*N*num_coeff, 1);
z_idx = zeros(M*N*num_coeff, 1);
one_over_num_frame = 1./(1+2*numFrame);
count = 1;

for j = 1:M
    for i = 1:N
        k = (j-1)*N + i;
        for f = -numFrame:numFrame
            
            jj = max(1, min(M, round(j + v(j, i)*f/numFrame)));
            ii = max(1, min(N, round(i + u(j, i)*f/numFrame)));
            kk = (jj - 1)*N + ii;
            x_idx(count) = k;
            y_idx(count) = kk;
            z_idx(count) = one_over_num_frame*0.9;
            count = count + 1;
            
            
%             jj = max(1, min(M, jj + 1));
%             ii = max(1, min(N, ii + 1));
% % 
% %             
%             kk = (jj - 1)*N + ii;
%             x_idx(count) = k;
%             y_idx(count) = kk;
%             z_idx(count) = one_over_num_frame*0.1;
%             count = count + 1;
%             
        end
    end
end
count = count - 1;



% K = sparse(x_idx(1:count), y_idx(1:count), z_idx(1:count), M*N, M*N);
K = sparse(x_idx(1:count), y_idx(1:count), 1./num_coeff, M*N, M*N);
% % K = sparse(y_idx(1:count), x_idx(1:count), z_idx, M*N, M*N);




%
%
% %
% [M, N] = size(u);
%
% numFrame = 75;%80;%100;
% one_over_num_frame = 1./numFrame;
%
% x_idx = zeros(M*N*(1+2*numFrame), 1);
% y_idx = zeros(M*N*(1+2*numFrame), 1);
% count = 1;
%
% for j = 1:M
%     for i = 1:N
%         k = (j-1)*N + i;
%         for f = -numFrame:numFrame
%             jj = round(j + v(j, i)*f*one_over_num_frame);
%             ii = round(i + u(j, i)*f*one_over_num_frame);
%             if 1 <= ii && ii <= N && 1 <= jj && jj <= M
%                 kk = (jj - 1)*N + ii;
%                 x_idx(count) = k;        % from
%                 y_idx(count) = kk;       % to
%                 count = count + 1;
%             end
%         end
%     end
% end
% count = count - 1;
%
% K = sparse(y_idx(1:count), x_idx(1:count), 1/(1+2*numFrame), M*N, M*N);





% 
% 
% 
% %
% 
% [M, N] = size(u);
% 
% numFrame = 70;%100;
% one_over_num_frame = 1./numFrame;
% 
% x_idx = zeros(M*N*(1+2*numFrame), 1);
% y_idx = zeros(M*N*(1+2*numFrame), 1);
% z_idx = zeros(M*N*(1+2*numFrame), 1);
% 
% 
% histogram = zeros(M*N, 1);
% count = 1;
% for j = 1:M
%     for i = 1:N
%         k = (j-1)*N + i;
%         for f = -numFrame:numFrame
%             jj = round(j + v(j, i)*f*one_over_num_frame);
%             ii = round(i + u(j, i)*f*one_over_num_frame);
%             if 0 < ii && ii <= N && 0 < jj && jj <= M
%                 kk = (jj - 1)*N + ii;
%                 x_idx(count) = k;        % from
%                 y_idx(count) = kk;       % to
%                 histogram(kk) = histogram(kk) + 1;
%                 count = count + 1;
%             end
%         end
%     end
% end
% count = count - 1;
% %
% K = sparse(y_idx(1:count), x_idx(1:count), 1/(1+2*numFrame), M*N, M*N);
% 
% 
% 
% % for i = 1 : count
% %     z_idx(i) = 1./histogram(y_idx(i));
% % end
% 
% % for learning
% % K = sparse(y_idx(1:count), x_idx(1:count), z_idx(1:count), M*N, M*N);
% 










