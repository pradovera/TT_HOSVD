function [W, Z, Ys] = find_range_hada_unfolding_TT(A, B, pos, Zold, Ys, r, p, dir)
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end

    d = A.order;
    if d ~= B.order
        error('The order of the two tensors must coincide')
    end
    if ~isequal(A.size, B.size)
        error('All the mode sizes of the two tensors must coincide')
    end
    
    if strcmpi(dir, 'right')
        Acore = permute(A.U{pos}, [3, 2, 1]);
        Bcore = permute(B.U{pos}, [3, 2, 1]);
        posright = (1:pos-1);
    else %if strcmpi(dir, 'left') 
        Acore = A.U{pos};
        Bcore = B.U{pos};
        posright = (d:-1:pos+1);
    end
    
    if ~isempty(Ys{posright(end)})
        numgen = r + p - size(Ys{posright(end)}, 3);
    else
        numgen = r + p;
    end
    if numgen > 0
        Ys = update_stored_range_hada_TT(A, B, Ys, numgen, posright, dir);
    end
    
    Y = multiply_hada_TT(Acore, Bcore, Zold, Ys{posright(end)}(:, :, 1:r+p), r+p);
    
    [W, ~] = qr(Y, 0);
    Z = left_project_hada_TT(Acore, Bcore, Zold, W);
    W = reshape(W, [size(Zold, 1), A.size(pos), size(W, 2)]);
end
