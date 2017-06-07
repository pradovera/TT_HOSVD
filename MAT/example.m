close all; clear; clc
run('../startup')
rng(42);

m = 1e4;
n = 1e4;
r = 20;
lambda = 3;
sigmas = lambda.^-(0:r-1);

A = gen_matrix_decay(m, n, sigmas);
% A = U * diag(sigmas) * V';

%%%%%%%%%%%%%%%%
%%%%% Truncation
%%%%%%%%%%%%%%%%
%% Gaussian
r_new = 15;
p = 3;
tic
[W, Z] = truncate_randomized(A, r_new, p);
toc
norm(A - W*Z', 'fro')

%% Rank-1
r_new = 15;
p = 3;
tic
[W, Z] = truncate_randomized_rank1(A, r_new, p);
toc
norm(A - W*Z', 'fro')

%% Uniform
r_new = 15;
p = 3;
tic
[W, Z] = truncate_randomized_uniform(A, r_new, p);
toc
norm(A - W*Z', 'fro')

%%%%%%%%%%%%%%
%%%%% Rounding
%%%%%%%%%%%%%%
%% Frobenius
epsilon = 1e-4;
q = 2;
tic
[W, Z] = round_randomized(A, epsilon, q);
toc
norm(A - W*Z', 'fro')

%% Frobenius, explicit norm
epsilon = 1e-4;
tic
[W, Z] = round_randomized_normA(A, epsilon);
toc
norm(A - W*Z', 'fro')

%% Spectral
epsilon = 1e-4;
q = 2;
tic
[W, Z] = round_randomized_2norm(A, epsilon, q);
toc
norm(A - W*Z')

%% Spectral, explicit norm
epsilon = 1e-4;
q = 2;
tic
[W, Z] = round_randomized_2norm_PM(A, epsilon, q);
toc
norm(A - W*Z')
