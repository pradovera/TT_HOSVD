function [A_star, Ys] = truncate_TT_randomized(A, r, p, dir, Ys)
    d = A.order;
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~exist('Ys', 'var')
        Ys = cell(1, d);
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end

    if numel(r) == 1
        r = [1, repmat(r, 1, d-1), 1];
    elseif numel(r) == d + 1
        r([1,end]) = 1;
    else
        error('Wrong format for the target TT-ranks. It must either be a scalar or a vector of length (A.order + 1)') 
    end
    A_star = A;
    
    if strcmpi(dir, 'right')
        for i = d:-1:2
            
            rEff = min(r(i), A_star.rank(i));
            [W, Z, Ys] = find_range_unfolding_TT(A_star, i, Ys, rEff, p, 'right');

            A_star.U{i} = reshape(W', [size(W, 2), A_star.size(i), A_star.rank(i + 1)]);
            A_star.U{i - 1} = tensorprod(A_star.U{i - 1}, Z, 3);
        end
    else %if strcmpi(dir, 'left') 
        for i = 1:d-1
            
            rEff = min(r(i+1), A_star.rank(i+1));
            [W, Z, Ys] = find_range_unfolding_TT(A_star, i, Ys, rEff, p, 'left');
            
            A_star.U{i} = reshape(W, [A_star.rank(i), A_star.size(i), size(W, 2)]);
            A_star.U{i + 1} = tensorprod(A_star.U{i + 1}, Z, 1);
        end
    end
end
