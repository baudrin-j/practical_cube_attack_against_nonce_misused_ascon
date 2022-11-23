/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : rounds_1_to_4.hpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the functions needed to compute THE FIRST 4 ROUNDS of ASCON
 *           The computation is not exhaustive, it only computes the part of the
 *           ANF which is needed.
*/

#ifndef ROUNDS_1_TO_4_HPP
#define ROUNDS_1_TO_4_HPP

#include <iostream>
#include <fstream>
#include <set>
#include <array>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <map>

// a monomial represented as a boolean vector of size 320
using monom = std::array<uint64_t, 5>;

// a coordinate is seen as a set of monomials
using coor = std::set<monom>;

// ASCON state made of 320 coordinates
using state = std::array<coor, 320>;

using uint = unsigned int;

// polynomial whose monomials can only be 1, bi*ci, bi, ci
using coefficient = std::array<uint64_t, 4>;
// polynomial whose variables are v_i and coefficients are polynomials in 1, bi*ci, bi, ci.
using poly_map = std::map<uint64_t, coefficient>;

const coor add_coor(const coor &c1, const coor &c2);
const coor mult_coor(const coor &c1, const coor &c2, const std::function<bool(const monom&)> &condition_mult);
const std::array<poly_map, 320> get_l4(const state &);

#endif /* ROUNDS_1_TO_4_HPP */
