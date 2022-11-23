/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : cube_sum.h
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : parallelized cube-sum computation
*/
#ifndef CUBE_SUM_H
#define CUBE_SUM_H

#include <iostream>
#include <vector>
#include <omp.h>

#include "permutation.h"

void cube_sum(uint64_t* partial_init, const uint &rounds, \
		const std::vector<uint> &cube_index, const bool &last_linlayer, const bool &cst);

#endif /* CUBE_SUM_H */
