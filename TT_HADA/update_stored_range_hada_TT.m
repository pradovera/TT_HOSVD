function Ys = update_stored_range_hada_TT(A, B, Ys, numgen, posright, dir)
%%% UPDATE Ys with new random samples.
% Ys is a cell array of size (A.order).
%
% If dir = 'left', each cell j contains an order 3 tensor of size
% A.rank(j) * B.rank(j) * k, with a sample on each mode-3 slice.
% Each sample is of the form $(A*B)_{\geq j}^\top(\omega_d \otimes ...
% \otimes \omega_j)$ for gaussian vectors \omega_j, ..., \omega_d
% (reshaped into a matrix of suitable size).
%
% If dir = 'right', the same but with size A.rank(j+1) * k and samples of
% the form $(A*B)_{\leq j}^{(j)}(\omega_j \otimes ... \otimes \omega_1)$.
%
% The construction is recursive.
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
    
    if ~isempty(Ys{posright(end)})
        numgenEff = numgen + size(Ys{posright(end)}, 3);
    else
        numgenEff = numgen;
    end
    
    for j = posright
        n = A.size(j);
        if strcmpi(dir, 'right')
            rankAl = A.rank(j + 1);
            rankBl = B.rank(j + 1);
            Acore = permute(A.U{j}, [3, 2, 1]);
            Bcore = permute(B.U{j}, [3, 2, 1]);
        else %if strcmpi(dir, 'left')
            rankAl = A.rank(j);
            rankBl = B.rank(j);
            Acore = A.U{j};
            Bcore = B.U{j};
        end
        if ~isempty(Ys{j})
            oldsize = size(Ys{j}, 3);
        else
            oldsize = 0;
        end
        Omega = randn(n, numgenEff - oldsize);
        if j == 1 || j == d
            Ynew = reshape(kr(Acore, Bcore) * Omega, [rankBl, rankAl, numgenEff - oldsize]);
        else
            if strcmpi(dir, 'right')
                Yr = Ys{j - 1}(:, :, oldsize+1:end);
            else %if strcmpi(dir, 'left')
                Yr = Ys{j + 1}(:, :, oldsize+1:end);
            end
            Ynew = zeros(rankBl, rankAl, numgenEff - oldsize);
            for l = 1:n
                Yaux1 = tensorprod(Yr, permute(Bcore(:, l, :), [1,3,2]), 1);
                Yaux2 = tensorprod(Yaux1, permute(Acore(:, l, :), [1,3,2]), 2);
                Ynew = Ynew + permute(bsxfun(@times, permute(Yaux2, [1,3,2]), Omega(l, :)), [1,3,2]);
            end
        end
        Ys{j} = cat(3, Ys{j}, Ynew);
    end
end
