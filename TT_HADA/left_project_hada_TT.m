function V = left_project_hada_TT(Xl, Yl, Vold, U)
    n = size(Xl, 2);
    if n ~= size(Yl, 2)
        error('The second mode size of the two tensors must coincide')
    end
    
    Xlperm = permute(Xl, [3, 1, 2]);
    Ylperm = permute(Yl, [3, 1, 2]);
    
    r = size(U, 2);
    U = reshape(U, [size(Vold, 1), n, r]);

    V1p = permute(tensorprod(U, Vold', 1), [1, 3, 2]);
    V = zeros(size(Yl, 3), size(Xl, 3), r);
    for j = 1:n
        V1t = reshape(V1p(:, :, j), [size(Yl, 1), size(Xl, 1), r]);
        V = V + tensorprod(tensorprod(V1t, Ylperm(:, :, j), 1), Xlperm(:, :, j), 2);
    end
    V = reshape(V, [size(Xl, 3) * size(Yl, 3), r])';
end
