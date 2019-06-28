%   TV-L1 optical flow
%
%   Author: Thomas Pock
%
%   If you use this file or package for your work, please refer to the
%   following papers:
% 
%   [1] Antonin Chambolle and Thomas Pock, A first-order primal-dual
%   algorithm with applications to imaging, Technical Report, 2010
%
%   License:
%     Copyright (C) 2011 Institute for Computer Graphics and Vision,
%                      Graz University of Technology
function rho = calc_rho(I1, I2, u, v, isValid)

[M N C] = size(I1);

idx = repmat([1:N], M,1);
idy = repmat([1:M]',1,N);
 
idxx = idx + u;
idyy = idy + v;
m = (idxx > N-1) | (idxx < 2) | (idyy > M-1) | (idyy < 2);
rho = abs(I1 - interp2(I2,idxx,idyy,'cubic'));
rho(m) = 0.0;

if nargin == 5
    rho(~isValid) = 0;
end

end