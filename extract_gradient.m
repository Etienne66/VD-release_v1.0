function [ dx dy ] = extract_gradient( src )

    [m_row, n_col] = size(src);
    dx = zeros(m_row, n_col);
    dy = zeros(m_row, n_col);
    
    for(i = 1 : m_row)        
        for(j = 1 : n_col)
           
            if(j ~= n_col)
                dx(i, j) = src(i, j) - src(i, j + 1);
            end
            
            if(i ~= m_row)
                dy(i, j) = src(i, j) - src(i + 1, j);            
        end        
    end
end

