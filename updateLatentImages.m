function [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = ...
    updateLatentImages(K, L, L_, B_Pyramid, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda, err)

n_imlist = length(L);

for iter = 1  : 5
    [L, L_, dual_spatial, dual_temporal_forward, dual_temporal_backward] = ...
        updateDuals2(K, L, L_, B_Pyramid, dual_spatial, dual_temporal_forward, dual_temporal_backward, flow_forward, flow_backward, num_temporal_neighbors, lambda);
    
    L = updatePrimals(L, K, B_Pyramid, num_temporal_neighbors, lambda, err);
    L_ = extrapolate(L, L_);
    
    for n = 1 : n_imlist
        out_name = sprintf('out1_%04d.png', n);
        imwrite(uint8(255*min(1, max(0, L{n}))), out_name);
    end
    
    fprintf('iter:%d\n', iter);
end

end



