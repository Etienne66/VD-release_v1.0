function [ dx_t_A ] = nablax_t_A( A )

    [m_row n_col] = size(A);
    dx_t_A = zeros(m_row, n_col);
    
    for(i = 1 : m_row)
        for(j = 1 : n_col)           
            if(j == 1)
                dx_t_A(i, j) = A(i, j);
            else if(j == n_col)
                dx_t_A(i, j) = - A(i, j - 1);
            else
                dx_t_A(i, j) = A(i, j) - A(i, j - 1);
            end
        end
    end


end

