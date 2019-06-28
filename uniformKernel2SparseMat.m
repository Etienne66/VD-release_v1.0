function [ K ] = uniformKernel2SparseMat( uniform_kernel, img_height, img_width )

    
img_dims = img_height*img_width;
kernel = flipud(fliplr(uniform_kernel));

[kernel_wdith kernel_height] = size(kernel);

valid_kernel = [];
valid_kernel_idx = [];

for(i = 1: 1 : kernel_height)    
    for(j = 1 : 1 : kernel_wdith)
        if( kernel(i, j) > 0 )
            valid_kernel = [valid_kernel kernel(i, j)];
%             valid_kernel_idx = [valid_kernel_idx ([i; j] - [ floor(kernel_height/2); floor(kernel_wdith/2) ] )];
valid_kernel_idx = [valid_kernel_idx ([i; j] - [ ceil(kernel_height/2); ceil(kernel_wdith/2) ] )];
        end
    end    
end

valid_kernel_dims = size(valid_kernel, 2);


% 
a = zeros(valid_kernel_dims*img_dims, 1);
b = zeros(valid_kernel_dims*img_dims, 1);
c = zeros(valid_kernel_dims*img_dims, 1);

for(i = ceil(kernel_height/2) : 1 : img_height - ceil(kernel_height/2)) 
    for(j = ceil(kernel_wdith/2) : 1 : img_width - ceil(kernel_wdith/2))

        for(k = 1 : 1 : valid_kernel_dims)

            idx = (k - 1)*img_dims + ((i - 1)*img_width + j);
            
            a(idx) = ( (i - 1)*img_width + j );            
            b(idx) = (i - 1 + valid_kernel_idx(1, k))*img_width + j + valid_kernel_idx(2, k);
            c(idx) = valid_kernel(k);
        end
    end
    
end

a(a == 0) = [];
b(b == 0) = [];
c(c == 0) = [];


K = sparse(a, b, c, img_dims, img_dims);

end

