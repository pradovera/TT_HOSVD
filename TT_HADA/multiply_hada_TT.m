function R = multiply_hada_TT(Xl, Yl, V, R0, r)
    n = size(Xl, 2);
    if n ~= size(Yl, 2)
        error('The second mode size of the two tensors must coincide')
    end
    
    R = zeros(size(V, 1), n, r);
    for j = 1:n
        Raux1 = tensorprod(R0, permute(Yl(:, j, :), [1,3,2]), 1);
        Raux2 = reshape(tensorprod(Raux1, permute(Xl(:, j, :), [1,3,2]), 2), [size(Yl, 1) * size(Xl, 1), r]);
        R(:, j, :) = V * Raux2;
    end
    R = reshape(R, [size(V, 1) * n, r]);
end

