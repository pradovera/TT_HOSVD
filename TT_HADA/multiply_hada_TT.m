function Y = multiply_hada_TT(Al, Bl, Z, Y0, r)
    n = size(Al, 2);
    if n ~= size(Bl, 2)
        error('The second mode size of the two tensors must coincide')
    end
    
    Y = zeros(size(Z, 1), n, r);
    for j = 1:n
        Yaux1 = tensorprod(Y0, permute(Bl(:, j, :), [1,3,2]), 1);
        Yaux2 = reshape(tensorprod(Yaux1, permute(Al(:, j, :), [1,3,2]), 2), [size(Bl, 1) * size(Al, 1), r]);
        Y(:, j, :) = Z * Yaux2;
    end
    Y = reshape(Y, [size(Z, 1) * n, r]);
end

