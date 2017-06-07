function [U, V, Rs] = find_range_unfolding_TT(X, pos, Rs, k, p, dir)
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end

    d = X.order;
    
    if strcmpi(dir, 'right')
        Al = unfold(X.U{pos}, 'right')';
        posright = (1:pos-1);
    else %if strcmpi(dir, 'left') 
        Al = unfold(X.U{pos}, 'left');
        posright = (d:-1:pos+1);
    end
    
    numgen = k + p - size(Rs{posright(end)}, 2);
    if numgen > 0
        Rs = update_stored_range_TT(X, Rs, numgen, posright, dir);
    end
    
    Y0 = Rs{posright(end)}(:, 1:k+p);
    Y = Al * Y0;

    [U, ~] = qr(Y, 0);
    V = U' * Al;
end
