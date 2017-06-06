function [ Wp, Zp ] = recompress(W, Z, ep)
    [U, S, V] = svd(Z, 'econ');
    
    summ = cumsum(diag(S).^2,'reverse');

    s = find(summ / summ(1) < ep^2, 1);
    if isempty(s)
        s = size(S, 1) + 1;
    end
    s = s - 1;
    
    Wp = W * V(:, 1:s);
    Zp = U(:, 1:s) * S(1:s, 1:s);
end

