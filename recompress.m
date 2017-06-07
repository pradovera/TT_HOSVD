function [ Wp, Zp ] = recompress(W, Z, epsilon)
    [U, S, V] = svd(Z, 'econ');
    
    summ = cumsum(diag(S).^2,'reverse');

    r = find(summ / summ(1) < epsilon^2, 1);
    if isempty(r)
        r = size(S, 1) + 1;
    end
    r = r - 1;
    
    Wp = W * V(:, 1:r);
    Zp = U(:, 1:r) * S(1:r, 1:r);
end

