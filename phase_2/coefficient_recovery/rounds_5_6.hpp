/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : rounds_5_6.hpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the functions needed to compute THE LAST 2 ROUNDS of ASCON
 *           The computation is not exhaustive, it only computes the part of the
 *           ANF which is needed.
*/

#ifndef ROUNDS_5_6_HPP
#define ROUNDS_5_6_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include "rounds_1_to_4.hpp"

using size_2_products = std::array<uint, 3>;
using generic_size_2_products = std::tuple<uint, uint>;
using trails = std::pair<size_2_products, size_2_products>;

const std::string coefficient_recovery(const uint &col, const std::array<poly_map, 320> &l4, const uint64_t &target);

#endif /* ROUNDS_5_6_HPP */
