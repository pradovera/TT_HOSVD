function [W, Z] = round_randomized(A, epsilon, q)
%%% Algorithm 3
    n = size(A, 2);
    epsilon2 = epsilon^2;
    if epsilon2 < 2 * eps
        warning('Effective precision smaller than machine epsilon');
        epsilon2 = 2 * eps;
    end
    
    Omega = randn(n, q);
    [W, ~] = qr(A * Omega, 0);
    Z = A' * W;
    
    normZ2 = cumsum(sum(Z.^2, 1))';
    
    j = 0;
    while (j == 0) || (((normZ2(j + q) - normZ2(j)) > epsilon2 * normZ2(j + q)) &&...
            ((j == 1) || (abs((normZ2(j + q) - normZ2(j)) / normZ2(j + q) -...
            (normZ2(j + q - 1) - normZ2(j - 1)) / normZ2(j + q - 1)) > 10 * eps)))
        j = j + 1;
        w = A * randn(n, 1);
        for ii = 1:2
            w = w - W * (w' * W)';
        end
        w = w / norm(w);
        W = [W, w];
        z = A' * w;
        Z = [Z, z];
        normZ2 = [normZ2; normZ2(end) + sum(z.^2)];
    end
end
