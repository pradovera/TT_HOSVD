function [Z, Rs] = truncate_hada_TT_randomized(X, Y, r, p, dir, recompress, Rs)
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~exist('recompress', 'var')
        recompress = false;
    end
    d = X.order;
    if d ~= Y.order
        error('The order of the two tensors must coincide')
    end
    if ~isequal(X.size, Y.size)
        error('All the mode sizes of the two tensors must coincide')
    end
    
    if ~exist('Rs', 'var')
        Rs = cell(1, d);
    end
    cores = cell(1, d);
    V = 1;

    if strcmpi(dir, 'right')
        for i = d:-1:2
            k = min(r(i), X.rank(i) * Y.rank(i));
            [U, V, Rs] = find_range_hada_unfolding_TT(X, Y, i, V, Rs, k, p, 'right');
            
            if recompress
                [U1, S1, V1] = svd(V, 'econ');
                k = trunc_singular(diag(S1), r);
                U = tensorprod(U, U1(:, 1:k)', 3);
                V = S1(1:k, 1:k) * V1(:, 1:k)';
            end

            cores{i} = permute(U, [3, 2, 1]);
        end
        i = 1;
        A = unfold(X.U{i}, 'left')';
        B = unfold(Y.U{i}, 'left')';
        VAB = V * kr(A, B);
        cores{i} = permute(reshape(VAB, [size(V, 1), X.size(i), 1]), [3, 2, 1]);
   elseif strcmpi(dir, 'left') 
        for i = 1:d-1
            k = min(r(i + 1), X.rank(i + 1) * Y.rank(i + 1));
            [U, V, Rs] = find_range_hada_unfolding_TT(X, Y, i, V, Rs, k, p, 'left');
            
            if recompress
                [U1, S1, V1] = svd(V, 'econ');
                k = trunc_singular(diag(S1), r);
                U = tensorprod(U, U1(:, 1:k)', 3);
                V = S1(1:k, 1:k) * V1(:, 1:k)';
            end

            cores{i} = U;
        end
        i = d;
        A = unfold(X.U{i}, 'right');
        B = unfold(Y.U{i}, 'right');
        VAB = V * kr(A, B);
        cores{i} = reshape(VAB, [size(V, 1), X.size(i), 1]);
    else
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
    Z = TTeMPS(cores);
end
