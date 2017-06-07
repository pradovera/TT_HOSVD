function [Ys] = update_stored_range_TT(A, Ys, numgen, posright, dir)
%%% UPDATE Ys with new random samples.
% Ys is a cell array of size (A.order).
%
% If dir = 'left', each cell j contains a matrix of size A.rank(j) * k,
% with a sample on each column. Each sample is of the form
% $A_{\geq j}^\top(\omega_d \otimes ... \otimes \omega_j)$ for gaussian
% vectors \omega_j, ..., \omega_d.
%
% If dir = 'right', the same but with size A.rank(j+1) * k and samples of
% the form $A_{\leq j}^{(j)}(\omega_j \otimes ... \otimes \omega_1)$.
%
% The construction is recursive.
    if ~exist('dir', 'var')
        dir = 'left';
    end
    if ~strcmpi(dir, 'right') && ~strcmpi(dir, 'left') 
        error('Unknown direction specified. Choose either LEFT or RIGHT') 
    end
    
    d = A.order;

    numgenEff = numgen + size(Ys{posright(end)}, 2);
    
    for j = posright
        if strcmpi(dir, 'right')
            rankl = A.rank(j + 1);
            Acore = permute(A.U{j}, [3, 2, 1]);
        else %if strcmpi(dir, 'left')
            rankl = A.rank(j);
            Acore = A.U{j};
        end
        oldsize = size(Ys{j}, 2);
        Omega = randn(numgenEff - oldsize, A.size(j));
        Y = tensorprod(Acore, Omega, 2);
        if j == 1 || j == d
            Ynew = squeeze(Y);
        else
            if strcmpi(dir, 'right')
                Yr = Ys{j - 1}(:, oldsize+1:end);
            else %if strcmpi(dir, 'left')
                Yr = Ys{j + 1}(:, oldsize+1:end);
            end
            Ynew = zeros(rankl, numgenEff - oldsize);
            for l = 1:numgenEff - oldsize
                Ynew(:, l) = permute(Y(:, l, :), [1, 3, 2]) * Yr(:, l);
            end
        end
        Ys{j} = [Ys{j}, Ynew];
    end
end
