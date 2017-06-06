function x = round_nonortho( x, tol, dir )
    %   ROUND algorithm from TTeMPS Toolbox.without orthogonalization step
    %   Michael Steinlechner, 2013-2016
    %   Davide Pradovera, 2017
    
    if ~exist('dir', 'var')
        dir = 'right';
    end

    sz = x.size;
    d = x.order;

    if strcmpi(dir, 'left')
        left_rank = 1;
        for i = 1:d-1
            [U,S,V] = svd( unfold( x.U{i}, 'left'), 'econ' );
            r = trunc_singular( diag(S), tol, true );
            U = U(:,1:r);
            V = V(:,1:r);
            S = S(1:r,1:r);
            x.U{i} = reshape( U, [left_rank, sz(i), r] );
            x.U{i+1} = tensorprod( x.U{i+1}, S*V', 1 );
            left_rank = r;
        end
    elseif strcmpi(dir, 'right') 
        right_rank = 1;
        for i = d:-1:2
            [U,S,V] = svd( unfold( x.U{i}, 'right'), 'econ' );
            r = trunc_singular( diag(S), tol, true );
            U = U(:,1:r);
            V = V(:,1:r);
            S = S(1:r,1:r);
            x.U{i} = reshape( V', [r, sz(i), right_rank] );
            x.U{i-1} = tensorprod( x.U{i-1}, (U*S).', 3 );
            right_rank = r;
        end
    else
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
end
