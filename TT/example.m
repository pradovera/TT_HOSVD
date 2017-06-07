close all; clear; clc
run('../startup')
rng(42);

d = 10;
n = 50;
r = 30;
lambda = 3;
sigmas = lambda.^-(0:r-1);

A = gen_TT_tensor_decay(d, n, sigmas);

%% Truncation
r_new = [1, repmat(15, 1, d-1), 1];
p = 5;
direction = 'left';
tic
[A_star, ~] = truncate_TT_randomized(A, r_new, p, direction);
toc
norm(A - A_star, 'fro')

%% Rounding
epsilon = 1e-6;
step = 3;
direction = 'left';
q = 2;
r_guess = [1, repmat(10, 1, d-1), 1];
tic
[A_star, ~] = round_TT_randomized(A, epsilon, step, direction, q, r_guess);
toc
norm(A - A_star, 'fro')

