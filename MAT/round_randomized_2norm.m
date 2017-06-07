function [W, Z] = round_randomized_2norm(A, epsilon, q)
%%% Algorithm 5
    normA2 = PM(A, 0)^2;

    n = size(A, 2);
    epsilon2 = epsilon^2;
    if epsilon2 < 2 * eps
        warning('Effective precision smaller than machine epsilon');
        epsilon2 = 2 * eps;
    end
    
    Y = A * randn(n, q);
    
    W = []; Z = [];
    while 200/pi*max(sum(Y.^2,1)) > epsilon2 * normA2
        k = randi(q);
        w = Y(:, k);
        w = w / norm(w);
        W = [W, w];
        z = A' * w;
        Z = [Z, z];
        Y(:, k) = A * randn(n, 1);
        for ii = 1:2
            Y(:, k) = Y(:, k) - W * (Y(:, k)' * W)';
        end
        for i = [1:k-1, k+1:q]
            for ii = 1:2
                Y(:, i) = Y(:, i) - w * (Y(:, i)' * w)';
            end
        end
    end
end
