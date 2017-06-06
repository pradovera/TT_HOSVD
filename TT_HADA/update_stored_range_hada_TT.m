function Rs = update_stored_range_hada_TT(X, Y, Rs, numgen, posright, dir)
%%% UPDATE Rs with new random samples.
% Rs is a cell array of size (X.order).
%
% If dir = 'left', each cell j contains an order 3 tensor of size
% X.rank(j) * Y.rank(j) * k, with a sample on each mode-3 slice.
% Each sample is of the form $(X*Y)_{\geq j}^\top(\omega_d \otimes ... \otimes
% \omega_j)$ for gaussian vectors \omega_j, ..., \omega_d (reshaped into a
% matrix of suitable size).
%
% If dir = 'right', the same but with size X.rank(j+1) * k and samples of
% the form $(X*Y)_{\leq j}^{(j)}(\omega_j \otimes ... \otimes \omega_1)$.
%
% The construction is recursive.
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
    
    if ~isempty(Rs{posright(end)})
        numgenEff = numgen + size(Rs{posright(end)}, 3);
    else
        numgenEff = numgen;
    end
    
    for j = posright
        n = X.size(j);
        if strcmpi(dir, 'right')
            rankXl = X.rank(j + 1);
            rankYl = Y.rank(j + 1);
            Xr = permute(X.U{j}, [3, 2, 1]);
            Yr = permute(Y.U{j}, [3, 2, 1]);
        elseif strcmpi(dir, 'left')
            rankXl = X.rank(j);
            rankYl = Y.rank(j);
            Xr = X.U{j};
            Yr = Y.U{j};
        end
        if ~isempty(Rs{j})
            oldsize = size(Rs{j}, 3);
        else
            oldsize = 0;
        end
        Omega = randn(n, numgenEff - oldsize);
        if j == 1 || j == d
            Rnew = reshape(kr(Xr, Yr) * Omega, [rankYl, rankXl, numgenEff - oldsize]);
        else
            if strcmpi(dir, 'right')
                Rr = Rs{j - 1}(:, :, oldsize+1:end);
            else%if strcmpi(dir, 'left')
                Rr = Rs{j + 1}(:, :, oldsize+1:end);
            end
            Rnew = zeros(rankYl, rankXl, numgenEff - oldsize);
            for l = 1:n
                Raux1 = tensorprod(Rr, permute(Yr(:, l, :), [1,3,2]), 1);
                Raux2 = tensorprod(Raux1, permute(Xr(:, l, :), [1,3,2]), 2);
                Rnew = Rnew + permute(bsxfun(@times, permute(Raux2, [1,3,2]), Omega(l, :)), [1,3,2]);
            end
        end
        Rs{j} = cat(3, Rs{j}, Rnew);
    end
end
