# Practical Cube Attack against Nonce-Misused Ascon - Supplementary Material

This repository contains supplementary material provided with [our paper](https://tosc.iacr.org/index.php/ToSC/article/view/xxxx) published in
Volume 2022, Issue 3 of [*IACR Transactions on Symmetric Cryptology*](https://tosc.iacr.org/).

It is organized as follows.

## `conditional_cubes_choices`

In this folder are stored 5 SageMath/Python scripts (we used SageMath Version 9.2).

- `ascon_functions.py` contains the linear layer and the s-box layer of ASCON. It is used by all the following files.
- `conditional_cubes_choices.py` contains justification for our choices of conditional cubes. The sets $S_i$ from our work are computed from the ANF of two rounds of ASCON.
- `conditional_cubes_counter_measures.py` explores some possible settings which could be used as counter measures. In those settings, no conditional cubes (or less effective ones) can be built in the way presented in our paper. This file follows the logic of file `conditional_cubes_choices.py`.



The two remaining files are provided to compare our method with two other works:

- `LDW17_conditional_cubes_choices.py` shows that in the previous work of Li, Dong & Wang [TOSC:LDW17], the way cubes are chosen for round-reduced initializations can actually be explained in the same way we present in our work.
- `CHK22_conditional_cubes_choices.py` splits the set of public variables into subsets. Those subsets give a more precise insight into the internal differences between our work and the work of Chang, Hong & Kang [EPRINT:CHK22].



## `phase_1`

This folder corresponds to the first phase of our attack: the recovery of most of the bits of vectors $a$ and $e$.

It contains the C++ files used to derive the results which underpin the assumptions introduced in our paper. A Makefile is also given. The subfolder `results` is voluntarily empty: it will contain the result files after computation.

- `phase_1_verification.cpp` is the file containing the main function.
- `cube_sum.cpp` provides a parallelized cube-sum function using OpenMP.
- `permutation.cpp` contains the permutation used in ASCON.
- `random.cpp` contains pseudo-random 64-bit word generation functions using the C++ standard library.



## `phase_2`

This folder corresponds to the second phase of our attack: the recovery of the remaining bits of $a$, when some of the bits of $a$ and all the bits of $e$ have already been recovered.

It contains C++ files, as well a SageMath script. The subfolder `results` is voluntarily empty: as above, it will contain the result files after computation. A Makefile is provided. By using `make phase_2` a file `phase_2.out` is compiled. Launching this program will launch the main function located in`coefficient_recovery/coefficient_recovery.cpp`.

The main loop of this function is repeated until there is no more bit left to recover or if the maximal number of tries is reached. It follows the high level steps listed below.

- First of all, a monomial of degree 32 whose coefficients depend on unknown variables $a_i$ is targeted. Then, we recover the polynomial expression of its coefficients after the sixth S-box layer.  This step uses functions from files `rounds_1_to_4.cpp` and `rounds_5_6.cpp`, which is  located in subfolder `coefficient_recovery`. The recovery updates an initial state and computes the necessary part of the ANF round after round. 
- Then, the corresponding cube-sum vector is computed. It uses function `cube_sum_given_cubes_given_a_e`  from file `values_recovery.cpp` and other auxiliary functions which are all located in files from folder `values_recovery`.
- Finally, the corresponding system is built and solved by calling the SageMath script `system_solving.py` at the root of this folder. If information can be recovered from this solving, then it is taken into account for the next loop.

/!\ NB : In order for the program to work properly three files have to be MODIFIED:

- the file `coefficient_recovery/coefficient_recovery.cpp` at line 200: change the shell, if necessary.
- the script `script.run` (auxiliary script used to call `system_solving.py`): some additional lines may be needed in order to correctly launch SageMath (for instance, we used the Conda package management system to handle Sage's dependencies).
- the script `system_solving.py` at line 24: `dir_results` needs to be set with the correct path leading to the subfolder `results` .



## `phase_3`

This folder corresponds to the last phase of our attack: the recovery of all bits of $b$ and $c$, when all bits of $a$ and $e$ have already been recovered.

It contains C++ files as well as a SageMath script. The subfolder `results` is voluntarily empty: it will contain the result files after computations.

It should be used as follows.

- First of all use the files in subfolder `coefficient_recovery` (a Makefile is provided inside the subfolder). The main function is present in file `coefficient_recovery.cpp`. It enables the recovery of coefficients of some degree-31 monomials after the sixth S-box layer  This file uses the main functions of files `rounds_1_to_4.cpp` and `rounds_5_6.cpp` to update an initial state and compute the necessary part of the ANF round after round. The results will be output in subfolder `results` as two different files.

  -  `parameters.txt` which will contain the pseudo-random values of $a$ and $e$, as well as the description of the targeted degree-31 monomials.
  - `polynomials_cube_x.txt` (where `x` is the index of the current targeted cube) in which all the 64 coefficients we will stored as polynomials in $b_i$ and $c_i$ bits.

These two files are used in the next steps.

Then, use files in subfolder `values_recovery`  (a Makefile is provided inside the subfolder). The main function is present in file `values_recovery.cpp`. It enables the recovery of the cube-sum vectors for each of the targeted degree-31 monomials provided thanks to file `parameters.txt`. Before computing the cube-sum, the values of $a$ and $e$ provided by `parameters.txt` are used, and some pseudo-random values for $b$ and $c$ are computed. Those newly-selected values, are stored in the subfolder `results`  in the new output file `cube_sum_vectors.txt`. Then, the newly-computed values of the cube-sum vectors are also stored in the same file. This file will be used in the next step.

NB : `values_recovery.cpp` uses the other C++ files present in the folder, namely `cube_sum.cpp` and `permutation.cpp` which are the same files as the ones present (and already presented) in folder `phase_1`.

Finally, the last step uses the SageMath script at the root of the folder `phase_3`: `system_solving.py`. This script uses `cube_sum_vectors.txt` as well as `polynomials_cube_x.txt` to build the corresponding system, to recover information from the system, and finally to verify if the information recovered is correct.



/!\ Phase 2 and 3 share a common framework, that is why files in both subfolders really look alike. However, we would like to emphasize that the differences between them are very important, as they enable the recovery of two disjoint sets of bits. We tried to emphasize as much as possible the differences between the two folders with comments.


## Authors
- [Jules Baudrin](https://who.paris.inria.fr/Jules.Baudrin/)
- [Anne Canteaut](https://www.rocq.inria.fr/secret/Anne.Canteaut/)
- [Léo Perrin](https://who.paris.inria.fr/Leo.Perrin/)



