function [ dy_t_A ] = nablay_t_A( A )

    [m_row n_col] = size(A);
    dy_t_A = zeros(m_row, n_col);
    
    for(i = 1 : m_row)
        for(j = 1 : n_col)           
            if(i == 1)
                dy_t_A(i, j) = A(i, j);
            else if(i == m_row)
                dy_t_A(i, j) = - A(i-1, j);
            else
                dy_t_A(i, j) = A(i, j) - A(i - 1, j);
            end
        end
    end


end

