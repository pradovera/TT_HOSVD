function [W, Z] = round_randomized_2norm_PM(A, epsilon, q)
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
        end
        normR = PM(A, Z);
    end
end
