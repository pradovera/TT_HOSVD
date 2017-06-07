function [ lambda ] = PM(A, Z, target)
%   approximate 2-norm of A-WZ'
    if ~exist('target', 'var')
        target = -1;
    end
    
    epsilon = 1e-3;
    
    q = randn(size(A, 2), 1);
    q = q / norm(q);
    lambda = 0;
    for i = 1:25
        y = A' * (A * q) - Z * (q' * Z)';
        lambdaOld = lambda;
        lambda = q' * y;
        if (target > 0 && lambda > target) || (lambda - lambdaOld < epsilon * lambda)
            lambda = lambda^.5;
            return;
        end
        q = y / norm(y);
    end
    lambda = lambda^.5;
end

