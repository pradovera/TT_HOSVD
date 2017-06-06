# TT_HOSVD

This repository contains the code used in the Master Thesis `Randomized low-rank approximation of matrices and tensors' by Davide Pradovera (davide.pradovera@epfl.ch).

## Main folder
The main folder contains the following files:
- LICENSE
- README.md
- PM.m (power method for general matrices)
- gen_TT_tensor_decay.m (algorithm to generate TT-tensors whose unfoldings have a prescribed decay)
- recompress.m (method to recompress a matrix in low-rank format)
- round_nonortho.m (rounding of a TT-tensor without orthogonalization step)
- truncate_nonortho.m (truncation of a TT-tensor without orthogonalization step)

## ./MAT folder
The MAT folder contains the files regarding the compression of matrices (Chapter 2 in the Thesis):
- round_randomized.m (rounding of a matrix using Algorithm 3 in the Thesis)
- round_randomized_2norm.m (rounding of a matrix using Algorithm 5 in the Thesis)
- round_randomized_2norm_PM.m (rounding of a matrix using Algorithm 4 in the Thesis)
- round_randomized_normA.m (rounding of a matrix using Algorithm 2 in the Thesis)
- truncate_randomized.m (truncation of a matrix using Algorithm 1 in the Thesis with Gaussian random vectors)
- truncate_randomized_rank1.m (truncation of a matrix using Algorithm 1 in the Thesis with rank-1 random vectors)
- truncate_randomized_uniform.m (truncation of a matrix using Algorithm 1 in the Thesis with uniform random vectors)

## ./MAT_HADA folder
The MAT_HADA folder contains the files regarding the compression of Hadamard products of matrices (Section 4.1 in the Thesis):
- round_hada_randomized.m (rounding of a Hadamard product of matrices using Algorithm 13 in the Thesis)
- truncate_hada_randomized.m (truncation of a Hadamard product of matrices using a modification of Algorithm 1 in the Thesis)

## ./TT folder
The TT folder contains the files regarding the compression of TT-tensors (Chapter 3 in the Thesis):
- find_range_unfolding_TT.m (single step within Algorithm 11 in the Thesis)
- round_TT_randomized.m (rounding of a TT-tensor using Algorithm 12 in the Thesis)
- truncate_TT_randomized.m (truncation of a TT-tensor using Algorithm 11 in the Thesis)
- update_stored_range_TT.m (algorithm for the extraction of new samples within Algorithms 11 and 12 in the Thesis)

## ./TT_HADA folder
The TT_HADA folder contains the files regarding the compression of Hadamard products of TT-tensors (Section 4.3 in the Thesis):
- find_range_hada_unfolding_TT.m (single step within Algorithm 14 in the Thesis)
- left_project_hada_TT.m (computation of line 10 within Algorithm 14 in the Thesis)
- multiply_hada_TT.m (computation of line 6 within Algorithm 14 in the Thesis)
- round_hada_TT_randomized.m (rounding of a Hadamard product of TT-tensors using Algorithm 15 in the Thesis)
- truncate_hada_TT_randomized.m (truncation of a Hadamard product of TT-tensors using Algorithm 14 in the Thesis)
- update_stored_range_hada_TT.m (algorithm for the extraction of new samples within Algorithms 14 and 15 in the Thesis)
