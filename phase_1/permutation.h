/*
 * Filename : permutation.h
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the necessary functions to apply ASCON permutation to a 320-bit
 * state. A state is represented as an array of 5 uint64_t.
*/
#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "stdint.h"
#include "stdio.h"

uint64_t rotr64(uint64_t, unsigned int);
void multi_p(uint64_t* x, unsigned int nb_rounds, bool cst, bool last_linlayer);

#endif /* PERMUTATION_H */
