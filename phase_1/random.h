/*
 * Filename : random.h
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Random uint64_t generation functions
*/
#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <stdint.h>

using monom = uint64_t;
using uint = unsigned int;

monom random_monom_weight(uint w);
monom random_monom();

#endif /* RANDOM_H */
