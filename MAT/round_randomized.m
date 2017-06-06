function [ W, Z ] = round_randomized(A, epsilon, q)
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
    
    normfro2 = cumsum(sum(Z.^2, 1))';
    
    j = 0;
    while (j == 0) || (((normfro2(j + q) - normfro2(j)) > epsilon2 * normfro2(j + q)) &&...
            ((j == 1) || (abs((normfro2(j + q) - normfro2(j)) / normfro2(j + q) -...
            (normfro2(j + q - 1) - normfro2(j - 1)) / normfro2(j + q - 1)) > 10 * eps)))
        j = j + 1;
        q0 = A * randn(n, 1);
        for ii = 1:2
            q0 = q0 - W * (q0' * W)';
        end
        q0 = q0 / norm(q0);
        W = [W, q0];
        s0 = A' * q0;
        Z = [Z, s0];
        normfro2 = [normfro2; normfro2(end) + sum(s0.^2)];
    end
end
