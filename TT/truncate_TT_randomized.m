function [Z, Rs] = truncate_TT_randomized(X, r, p, dir, Rs)
    d = X.order;
    Z = X;
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~exist('Rs', 'var')
        Rs = cell(1, d);
    end


    if strcmpi(dir, 'right')
        for i = d:-1:2
            
            k = min(r(i), Z.rank(i));
            [U, V, Rs] = find_range_unfolding_TT(Z, i, Rs, k, p, 'right');

            Z.U{i} = reshape(U', [size(U, 2), Z.size(i), Z.rank(i + 1)]);
            Z.U{i - 1} = tensorprod(Z.U{i - 1}, V, 3);
        end
    elseif strcmpi(dir, 'left') 
        for i = 1:d-1
            
            k = min(r(i+1), Z.rank(i+1));
            [U, V, Rs] = find_range_unfolding_TT(Z, i, Rs, k, p, 'left');
            
            Z.U{i} = reshape(U, [Z.rank(i), Z.size(i), size(U, 2)]);
            Z.U{i + 1} = tensorprod(Z.U{i + 1}, V, 1);
        end
    else
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
end
