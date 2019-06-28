%   generate the sparse matrix of gradient operator(nabla)

%   Authors: Thomas Pock and Yunjin Chen
%   License:
%     Copyright (C) 2011 Institute for Computer Graphics and Vision,
%                      Graz University of Technology
function [nabla] = make_nabla(M,N)
%%% o----->x(M) %%%%%%%%%%%%%%%%%%%%%%%%%
%%% |
%%% |
%%% ¡ýy(N)      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kxf{matrix for gradient_x}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;

for y=1:M-1
  for x=1:N
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M;
    val(cnt) = -1;
    cnt = cnt+1;
    
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M+1;
    val(cnt) = 1;
    cnt = cnt+1;
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);

Kxf = sparse(row,col,val,M*N,M*N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kyf(matrix for gradient_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
row  = zeros(1,M*N*2);
col = zeros(1,M*N*2);
val  = zeros(1,M*N*2);

cnt = 1;

for y=1:M
  for x=1:N-1
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x-1)*M;
    val(cnt) = -1;
    cnt = cnt+1;
    
    row(cnt) = y+(x-1)*M;
    col(cnt) = y+(x)*M;
    val(cnt) = 1;
    cnt = cnt+1;
  end
end
row = row(1:cnt-1);
col = col(1:cnt-1);
val = val(1:cnt-1);

Kyf = sparse(row,col,val,M*N,M*N);

nabla = [Kxf;Kyf];