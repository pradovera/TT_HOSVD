function [Rs] = update_stored_range_TT(X, Rs, numgen, posright, dir)
%%% UPDATE Rs with new random samples.
% Rs is a cell array of size (X.order).
%
% If dir = 'left', each cell j contains a matrix of size X.rank(j) * k,
% with a sample on each column. Each sample is of the form
% $X_{\geq j}^\top(\omega_d \otimes ... \otimes \omega_j)$ for gaussian
% vectors \omega_j, ..., \omega_d.
%
% If dir = 'right', the same but with size X.rank(j+1) * k and samples of
% the form $X_{\leq j}^{(j)}(\omega_j \otimes ... \otimes \omega_1)$.
%
% The construction is recursive.
    if ~exist('dir', 'var')
        dir = 'left';
    end
    d = X.order;

    numgenEff = numgen + size(Rs{posright(end)}, 2);
    
    for j = posright
        if strcmpi(dir, 'right')
            rankl = X.rank(j + 1);
            Ar = permute(X.U{j}, [3, 2, 1]);
        elseif strcmpi(dir, 'left')
            rankl = X.rank(j);
            Ar = X.U{j};
        end
        oldsize = size(Rs{j}, 2);
        Omega = randn(numgenEff - oldsize, X.size(j));
        R = tensorprod(Ar, Omega, 2);
        if j == 1 || j == d
            Rnew = squeeze(R);
        else
            if strcmpi(dir, 'right')
                Rr = Rs{j - 1}(:, oldsize+1:end);
            else%if strcmpi(dir, 'left')
                Rr = Rs{j + 1}(:, oldsize+1:end);
            end
            Rnew = zeros(rankl, numgenEff - oldsize);
            for l = 1:numgenEff - oldsize
                Rnew(:, l) = permute(R(:, l, :), [1, 3, 2]) * Rr(:, l);
            end
        end
        Rs{j} = [Rs{j}, Rnew];
    end
end
