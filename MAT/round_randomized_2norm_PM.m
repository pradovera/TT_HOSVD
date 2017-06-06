function [ W, Z ] = round_randomized_2norm_PM(A, epsilon, q)
%%% Algorithm 4
    normA = PM(A, 0);

    n = size(A, 2);
    if epsilon < 2 * eps^.5
        warning('Effective precision smaller than machine epsilon');
        epsilon = 2 * eps^.5;
    end
    
    W = []; Z = [];
    j = 0;
    normR = normA;
    while normR > epsilon * normA
        for i = 1:q
            j = j + 1;
            q0 = A * randn(n, 1);
            if ~isempty(W)
                for ii = 1:2
                    q0 = q0 - W * (q0' * W)';
                end
            end
            q0 = q0 / norm(q0);
            W = [W, q0];
            z0 = A' * q0;
            Z = [Z, z0];
        end
        normR = PM(A, Z);
    end
end
