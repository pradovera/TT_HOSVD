function A_star = round_nonortho(A, tol, dir)
%   ROUND algorithm from TTeMPS Toolbox.without orthogonalization step
%   Michael Steinlechner, 2013-2016
%   Davide Pradovera, 2017
    if ~exist('dir', 'var')
        dir = 'right';
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end

    d = A.order;
    A_star = A;
    
    if strcmpi(dir, 'right')
        right_rank = 1;
        for i = d:-1:2
            [U,S,V] = svd(unfold(A_star.U{i}, 'right'), 'econ');
            r = trunc_singular(diag(S), tol, true);
            U = U(:,1:r);
            V = V(:,1:r);
            S = S(1:r,1:r);
            A_star.U{i} = reshape(V', [r, A_star.size(i), right_rank]);
            A_star.U{i-1} = tensorprod(A_star.U{i-1}, (U*S).', 3);
            right_rank = r;
        end
    else %if strcmpi(dir, 'left') 
        left_rank = 1;
        for i = 1:d-1
            [U,S,V] = svd(unfold(A_star.U{i}, 'left'), 'econ');
            r = trunc_singular(diag(S), tol, true);
            U = U(:,1:r);
            V = V(:,1:r);
            S = S(1:r,1:r);
            A_star.U{i} = reshape(U, [left_rank, A_star.size(i), r]);
            A_star.U{i+1} = tensorprod(A_star.U{i+1}, S*V', 1);
            left_rank = r;
        end
    end
end
