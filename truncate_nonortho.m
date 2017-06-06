function x = truncate_nonortho( x, r, dir )
    %   TRUNCATE algorithm from TTeMPS Toolbox.without orthogonalization step
    %   Michael Steinlechner, 2013-2016
    %   Davide Pradovera, 2017
    
    if ~exist('dir', 'var')
        dir = 'right';
    end

    if strcmpi(dir, 'left')
        for i = 1:x.order-1
            [U,S,V] = svd( unfold( x.U{i}, 'left'), 'econ' );
            s = min( r(i + 1), length(S));
            U = U(:,1:s);
            V = V(:,1:s);
            S = S(1:s,1:s);
            x.U{i} = reshape( U, [x.rank(i), x.size(i), s] );
            x.U{i+1} = tensorprod( x.U{i+1}, S*V', 1 );
        end
    elseif strcmpi(dir, 'right') 
        for i = x.order:-1:2
            [U,S,V] = svd( unfold( x.U{i}, 'right'), 'econ' );
            s = min( r(i), length(S));
            U = U(:,1:s);
            V = V(:,1:s);
            S = S(1:s,1:s);
            x.U{i} = reshape( V', [s, x.size(i), x.rank(i+1)] );
            x.U{i-1} = tensorprod( x.U{i-1}, (U*S)', 3 );
        end
    else
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
end
