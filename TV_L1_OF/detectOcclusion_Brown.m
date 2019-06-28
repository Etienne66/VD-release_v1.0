function [ occ ] = detectOcclusion_Brown(u, v)

[M N] = size(u);

occCount = zeros(M, N);
occ = zeros(M, N);

for i = 1 : M
    for j = 1 : N
        x = max(1, min(N, round(j + u(i, j))));
        y = max(1, min(M, round(i + v(i, j))));
        occCount(y,x) = occCount(y,x) + 1;
    end
end

for i = 1 : M
    for j = 1 : N
        org_x = round(j + u(i, j));
        org_y = round(i + v(i, j));
        x = max(1, min(N, round(j + u(i, j))));
        y = max(1, min(M, round(i + v(i, j))));
        
        if org_x < 1 || org_x > N || org_y < 1 || org_y > M
            occ(i, j) = 0.01;
        elseif occCount(y,x) > 2
            occ(i, j) = 0.01;
        elseif occCount(y,x) == 2
            occ(i, j) = 0.5;
        else
            occ(i, j) = 1;
        end
    end
end


end
