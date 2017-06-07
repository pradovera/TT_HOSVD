function Z = left_project_hada_TT(Al, Bl, Zold, W)
    n = size(Al, 2);
    if n ~= size(Bl, 2)
        error('The second mode size of the two tensors must coincide')
    end
    
    Xlperm = permute(Al, [3, 1, 2]);
    Ylperm = permute(Bl, [3, 1, 2]);
    
    r = size(W, 2);
    W = reshape(W, [size(Zold, 1), n, r]);

    Z1p = permute(tensorprod(W, Zold', 1), [1, 3, 2]);
    Z = zeros(size(Bl, 3), size(Al, 3), r);
    for j = 1:n
        Z1t = reshape(Z1p(:, :, j), [size(Bl, 1), size(Al, 1), r]);
        Z = Z + tensorprod(tensorprod(Z1t, Ylperm(:, :, j), 1), Xlperm(:, :, j), 2);
    end
    Z = reshape(Z, [size(Al, 3) * size(Bl, 3), r])';
end
