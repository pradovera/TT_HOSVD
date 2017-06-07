function [W, Z] = round_randomized_normA(A, epsilon, normA2)
%%% Algorithm 2
    if ~exist('normA2', 'var')
        normA2 = norm(A, 'fro')^2;
    end
    
    n = size(A, 2);
    epsilon2 = epsilon^2;
    if epsilon2 < 2 * eps
        warning('Effective precision smaller than machine epsilon');
        epsilon2 = 2 * eps;
    end
    
    normZ2 = 0;
    W = []; Z = [];
    
    while normA2 - normZ2 > epsilon2 * normA2
        w = A * randn(n, 1);
        if ~isempty(W)
            for ii = 1:2
                w = w - W * (w' * W)';
            end
        end
        w = w / norm(w);
        W = [W, w];
        z = A' * w;
        Z = [Z, z];
        normZ2 = normZ2 + norm(z)^2;
    end
end
