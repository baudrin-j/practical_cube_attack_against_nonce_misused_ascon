/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : values_recovery.h
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Recovery of the cube-sum vector corresponding to a "parameter" file
 *           which was output by coefficient_recovery
*/
 
#ifndef VALUES_RECOVERY_H
#define VALUES_RECOVERY_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <omp.h>

#include "cube_sum.h"

void cube_sum_given_cubes_given_a_e(const std::string &inputfilename, const std::string &outputfilename);
uint64_t random_monom();

#endif /* VALUES_RECOVERY_H */
