function [L_] = extrapolate(L, L_)

num_imgs = length(L);
for n = 1 : num_imgs
    L_{n} = 2*L{n} - L_{n};
end



end


