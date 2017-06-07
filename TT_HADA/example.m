close all; clear; clc
run('../startup')
rng(42);

d = 10;
n = 50;
rA = 30;
rB = 20;
lambdaA = 3;
lambdaB = 5;
sigmasA = lambdaA.^-(0:rA-1);
sigmasB = lambdaB.^-(0:rB-1);

A = gen_TT_tensor_decay(d, n, sigmasA);
B = gen_TT_tensor_decay(d, n, sigmasB);

C = hadamard(A, B);

%% Truncation
r_new = [1, repmat(50, 1, d-1), 1];
p = 5;
direction = 'left';
tic
[C_star, ~] = truncate_hada_TT_randomized(A, B, r_new, p, direction);
toc
norm(C - C_star, 'fro')

%% Rounding
epsilon = 1e-6;
step = 3;
direction = 'left';
q = 2;
r_guess = [1, repmat(20, 1, d-1), 1];
tic
[C_star, ~] = round_hada_TT_randomized(A, B, epsilon, step, direction, q, r_guess);
toc
norm(C - C_star, 'fro')
