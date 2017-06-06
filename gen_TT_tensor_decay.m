function [ Z ] = gen_TT_tensor_decay( d, n, s )
    [Q, ~] = qr(randn(n, length(s)), 0);
    cores = cell(1, d);
    for j = 1:length(s)
        a = Q(:, j);
        for i = 1:d
            cores{i} = reshape(a, [1, n, 1]);
        end
        if j == 1
            Z = s(j) * TTeMPS(cores);
        else
            Z = Z + s(j) * TTeMPS(cores);
        end
    end
end

