function [C, Ys] = truncate_hada_TT_randomized(A, B, r, p, dir, Ys)
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
    
    if ~exist('Rs', 'var')
        Ys = cell(1, d);
    end
    cores = cell(1, d);
    Z = 1;

    if strcmpi(dir, 'right')
        for i = d:-1:2
            rEff = min(r(i), A.rank(i) * B.rank(i));
            [W, Z, Ys] = find_range_hada_unfolding_TT(A, B, i, Z, Ys, rEff, p, 'right');
            
            cores{i} = permute(W, [3, 2, 1]);
        end
        i = 1;
        Acore = unfold(A.U{i}, 'left')';
        Bcore = unfold(B.U{i}, 'left')';
        ZAB = Z * kr(Acore, Bcore);
        cores{i} = permute(reshape(ZAB, [size(Z, 1), A.size(i), 1]), [3, 2, 1]);
   else %if strcmpi(dir, 'left') 
        for i = 1:d-1
            rEff = min(r(i + 1), A.rank(i + 1) * B.rank(i + 1));
            [W, Z, Ys] = find_range_hada_unfolding_TT(A, B, i, Z, Ys, rEff, p, 'left');
            
            cores{i} = W;
        end
        i = d;
        Acore = unfold(A.U{i}, 'right');
        Bcore = unfold(B.U{i}, 'right');
        ZAB = Z * kr(Acore, Bcore);
        cores{i} = reshape(ZAB, [size(Z, 1), A.size(i), 1]);
    end
    C = TTeMPS(cores);
end
