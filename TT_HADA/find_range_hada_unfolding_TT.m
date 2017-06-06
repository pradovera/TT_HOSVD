function [U, V, Rs] = find_range_hada_unfolding_TT(X, Y, pos, Vold, Rs, k, p, dir)
    if ~exist('dir', 'var')
        dir = 'left';
    end
    d = X.order;
    if d ~= Y.order
        error('The order of the two tensors must coincide')
    end
    if ~isequal(X.size, Y.size)
        error('All the mode sizes of the two tensors must coincide')
    end
    
    if strcmpi(dir, 'right')
        Xl = permute(X.U{pos}, [3, 2, 1]);
        Yl = permute(Y.U{pos}, [3, 2, 1]);
        posright = (1:pos-1);
    elseif strcmpi(dir, 'left') 
        Xl = X.U{pos};
        Yl = Y.U{pos};
        posright = (d:-1:pos+1);
    else
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
    
    if ~isempty(Rs{posright(end)})
        numgen = k + p - size(Rs{posright(end)}, 3);
    else
        numgen = k + p;
    end
    if numgen > 0
        Rs = update_stored_range_hada_TT(X, Y, Rs, numgen, posright, dir);
    end
    
    R = multiply_hada_TT(Xl, Yl, Vold, Rs{posright(end)}(:, :, 1:k+p), k+p);
    
    [U, ~] = qr(R, 0);
    V = left_project_hada_TT(Xl, Yl, Vold, U);
    U = reshape(U, [size(Vold, 1), X.size(pos), size(U, 2)]);
end
