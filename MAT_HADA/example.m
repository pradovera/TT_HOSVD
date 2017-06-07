close all; clear; clc
run('../startup')
rng(42);

m = 1e4;
n = 1e4;
rA = 20;
rB = 15;
lambdaA = 3;
lambdaB = 5;
sigmasA = lambdaA.^-(0:rA-1);
sigmasB = lambdaB.^-(0:rB-1);

[A, UA, ~, VA] = gen_matrix_decay(m, n, sigmasA);
% A = UA * diag(sigmasA) * VA';
WA = bsxfun(@times, UA, sigmasA);
ZA = VA;


[B, UB, ~, VB] = gen_matrix_decay(m, n, sigmasB);
% B = UB * diag(sigmasB) * VB';
WB = bsxfun(@times, UB, sigmasB);
ZB = VB;

%% Truncation
r_new = 30;
p = 3;
tic
[W, Z] = truncate_hada_randomized(WA, ZA, WB, ZB, r_new, p);
toc
norm(A.*B - W*Z', 'fro')

%% Rounding
epsilon = 1e-4;
q = 1;
tic
[W, Z] = round_hada_randomized(WA, ZA, WB, ZB, epsilon, q);
toc
norm(A.*B - W*Z', 'fro')

%% Rounding with recompression
epsilon = 1e-4;
q = 1;
tic
[W0, Z0] = round_hada_randomized(WA, ZA, WB, ZB, epsilon, q);
[W, Z] = recompress(W0, Z0, epsilon);
toc
norm(A.*B - W*Z', 'fro')

