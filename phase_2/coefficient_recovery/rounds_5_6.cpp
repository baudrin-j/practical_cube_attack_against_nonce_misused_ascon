/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : rounds_5_6.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the functions needed to compute THE LAST 2 ROUNDS of ASCON
 *           The computation is not exhaustive, it only computes the part of the
 *           ANF which is needed.
*/

#include "rounds_5_6.hpp"

using namespace std;
using namespace chrono;


/*
 * Converts a coefficient to a readable string.
 */
string convert_monom_to_txt(const coor &c) {
	string s;
	bool plus = false;

	for(const auto &m: c) {
		if(!plus)
			plus = true;
		else
			s += " + ";

		bool times = false;
		for(uint j = 0; j < 64; j++) {
			if((m[1] >> (63 - j)) & 1) {
				if(!times)
					times = true;
				else
					s += "*";

				s += "a" + to_string(j);
			}
		}

		if(!times) // Handles the case of the cst coefficient
			s = "1";
	}
	if(s.empty())
		s = "0";
	return s;
}


/*
 * Computes a partial multiplication between two poly_map and returns a poly_map.
 * It is expected that c1 and c2 are two poly_maps with terms of degree 8, as output by get_l4.
 * It only returns the terms of degree 16 that appears in the genuine product.
 *
 * This function corresponds to the computation of the interesting terms during S5.
 */
const poly_map multiply_maps_S5(const poly_map &c1, const poly_map &c2) {
	poly_map prod; // Output product
	const function<bool(const monom &)> f_always_true = [](const monom &m) { return true; };

	// Double for loop to compute the product, restricted by no condition
	for(const auto &[monom1, coeff1]: c1) {
		for(const auto &[monom2, coeff2]: c2) {
			const uint64_t tmp_monom = (monom1 | monom2); // Multiplication of two monomials is an OR
			if(((uint) __builtin_popcountll(tmp_monom)) == 16) {
				prod[tmp_monom] = add_coor(prod[tmp_monom], mult_coor(coeff1, coeff2, f_always_true));
			}
		}
	}
	return prod;
}


/*
 * Computes a partial multiplication between two poly_map and returns a coefficient.
 * It is expected that c1 and c2 are two poly_maps with terms of degree 16, as output by multiply_maps_S5.
 * It only returns the coefficient (that appears in the genuine product c1*c2) corresponding to the target monomial given as input.
 * It is expected that target is of degree (hamming weight) 32.
 *
 * This function corresponds to the computation of a coefficient of a monomial of degree 32 after S6.
 */
coor multiply_maps_S6(const poly_map &c1, const poly_map &c2, const uint64_t &target) {
	coor prod; // Output coefficient
	const function<bool(const monom &)> f_always_true = [](const monom &m) { return true; };

	// Select the smallest list to be browsed
	const poly_map * first = &c1;
	const poly_map * second = &c2;
	if(c2.size() < c1.size()) {
		first = &c2;
		second = &c1;
	}

	for(const auto &[monom1, coeff1]: (*first)) { // Loop over the smallest list
		if(!coeff1.empty()) { // If monom1 actually appears
			const uint64_t complement = ((~monom1) & target);

			if((*second).contains(complement) && !(*second).at(complement).empty()) { // Look for the complementary monom in the second list
				prod = add_coor(prod, mult_coor(coeff1, (*second).at(complement), f_always_true ));
			}
		}
	}
	return prod;
}


/*
 * Computes the coefficient of a monomial of degree 32 in a single coordinate c_{0,x} after S6.
 * The coefficient is output as a string.
 *
 * - col is the index of the coordinate in which we are looking for. (0 <= col <= 63)
 * - target is the monomial we are targeting. (target is a word of size 64 and hamming weight 32)
 * - l4 is the state after l4, initialized with only the necessary variables.
 *  It is expected that l4 has been initialized through get_l4 first.
 *
 */
const string coefficient_recovery(const uint &col, const array<poly_map, 320> &l4, const uint64_t &target) {
	const auto start_s5 = high_resolution_clock::now();
	cout << "S5-L5..." << endl;

	// List of products of size 2 appearing in at least one trail of size 4
	const vector<size_2_products> list_products = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {0, 3, 4}, \
	{3, 1, 2}, {3, 1, 3}, {3, 2, 3}, {23, 0, 1}, {23, 1, 4}, {25, 1, 2}, {25, 1, 3}, {25, 2, 3}, \
	{36, 0, 1}, {36, 1, 2}, {36, 1, 4}, {45, 0, 1}, {45, 1, 2}, {45, 1, 4}, {57, 0, 1}, {57, 1, 4}, {58, 3, 4}, {63, 3, 4}};

	// Table mapping a product to its actual polynomial
	map<const size_2_products, poly_map> same_col_products;

	// STEP 1 : for each product of size 2, computes the product and store it in the table
#pragma omp parallel for default(none) shared(l4, list_products, same_col_products, std::cout, col)
	for(const auto &cur_prod : list_products) {
		const auto &[x, y1, y2] = cur_prod;
		const poly_map &c1 = l4[y1 * 64 + ((x + col) % 64)];
		const poly_map &c2 = l4[y2 * 64 + ((x + col) % 64)];

		cout << "Prod [" + to_string(x) + ", " + to_string(y1) + ", " + to_string(y2) +  "] - Nb checks:" + to_string((c1.size() * c2.size()) / 1000000) + "M\n";

		same_col_products[cur_prod] = multiply_maps_S5(c1, c2);
	}

	const auto stop_s5 = high_resolution_clock::now();
	const auto duration_s5 = duration_cast<seconds>(stop_s5 - start_s5);
	cout << "S5-L5 done in " + to_string(duration_s5.count()) + "secs.\nS6...";

	// List of trails of size 4 leading to coordinate c_{0,0} through 1.5 round
	const vector<trails> list_trails = {{{25, 2, 3}, {63, 3, 4}}, {{25, 2, 3}, {58, 3, 4}}, {{0, 3, 4}, {25, 2, 3}}, {{3, 2, 3}, {63, 3, 4}}, {{58, 3, 4}, {3, 2, 3}}, {{0, 3, 4}, {3, 2, 3}}, {{0, 2, 3}, {63, 3, 4}}, {{0, 2, 3}, {58, 3, 4}}, {{57, 1, 4}, {25, 2, 3}}, {{57, 1, 4}, {3, 2, 3}}, {{0, 2, 3}, {57, 1, 4}}, {{25, 2, 3}, {45, 1, 4}}, {{25, 2, 3}, {45, 1, 2}}, {{3, 2, 3}, {45, 1, 4}}, {{3, 2, 3}, {45, 1, 2}}, {{0, 2, 3}, {45, 1, 4}}, {{0, 2, 3}, {45, 1, 2}}, {{25, 2, 3}, {36, 1, 4}}, {{25, 2, 3}, {36, 1, 2}}, {{3, 2, 3}, {36, 1, 4}}, {{3, 2, 3}, {36, 1, 2}}, {{0, 2, 3}, {36, 1, 4}}, {{0, 2, 3}, {36, 1, 2}}, {{25, 1, 3}, {63, 3, 4}}, {{25, 1, 3}, {58, 3, 4}}, {{0, 3, 4}, {25, 1, 3}}, {{25, 1, 2}, {63, 3, 4}}, {{25, 1, 2}, {58, 3, 4}}, {{0, 3, 4}, {25, 1, 2}}, {{25, 1, 3}, {57, 1, 4}}, {{25, 1, 2}, {57, 1, 4}}, {{25, 1, 3}, {45, 1, 4}}, {{25, 1, 3}, {45, 1, 2}}, {{25, 1, 2}, {45, 1, 4}}, {{25, 1, 2}, {45, 1, 2}}, {{25, 1, 3}, {36, 1, 4}}, {{25, 1, 3}, {36, 1, 2}}, {{25, 1, 2}, {36, 1, 4}}, {{25, 1, 2}, {36, 1, 2}}, {{25, 2, 3}, {23, 1, 4}}, {{3, 2, 3}, {23, 1, 4}}, {{0, 2, 3}, {23, 1, 4}}, {{25, 1, 3}, {23, 1, 4}}, {{25, 1, 2}, {23, 1, 4}}, {{3, 1, 3}, {63, 3, 4}}, {{58, 3, 4}, {3, 1, 3}}, {{0, 3, 4}, {3, 1, 3}}, {{3, 1, 2}, {63, 3, 4}}, {{58, 3, 4}, {3, 1, 2}}, {{0, 3, 4}, {3, 1, 2}}, {{57, 1, 4}, {3, 1, 3}}, {{57, 1, 4}, {3, 1, 2}}, {{3, 1, 3}, {45, 1, 4}}, {{3, 1, 3}, {45, 1, 2}}, {{3, 1, 2}, {45, 1, 4}}, {{3, 1, 2}, {45, 1, 2}}, {{3, 1, 3}, {36, 1, 4}}, {{3, 1, 3}, {36, 1, 2}}, {{3, 1, 2}, {36, 1, 4}}, {{3, 1, 2}, {36, 1, 2}}, {{3, 1, 3}, {23, 1, 4}}, {{3, 1, 2}, {23, 1, 4}}, {{0, 1, 3}, {63, 3, 4}}, {{0, 1, 3}, {58, 3, 4}}, {{0, 1, 2}, {63, 3, 4}}, {{0, 1, 2}, {58, 3, 4}}, {{0, 1, 2}, {0, 3, 4}}, {{0, 1, 2}, {25, 2, 3}}, {{0, 1, 2}, {3, 2, 3}}, {{0, 1, 3}, {57, 1, 4}}, {{0, 1, 2}, {57, 1, 4}}, {{0, 1, 3}, {45, 1, 4}}, {{0, 1, 3}, {45, 1, 2}}, {{0, 1, 2}, {45, 1, 4}}, {{0, 1, 2}, {45, 1, 2}}, {{0, 1, 3}, {36, 1, 4}}, {{0, 1, 3}, {36, 1, 2}}, {{0, 1, 2}, {36, 1, 4}}, {{0, 1, 2}, {36, 1, 2}}, {{0, 1, 2}, {25, 1, 3}}, {{0, 1, 2}, {25, 1, 2}}, {{0, 1, 3}, {23, 1, 4}}, {{0, 1, 2}, {23, 1, 4}}, {{0, 1, 2}, {3, 1, 3}}, {{0, 1, 2}, {3, 1, 2}}, {{57, 0, 1}, {25, 2, 3}}, {{57, 0, 1}, {3, 2, 3}}, {{0, 2, 3}, {57, 0, 1}}, {{57, 0, 1}, {25, 1, 3}}, {{57, 0, 1}, {25, 1, 2}}, {{57, 0, 1}, {3, 1, 3}}, {{57, 0, 1}, {3, 1, 2}}, {{0, 1, 3}, {57, 0, 1}}, {{0, 1, 2}, {57, 0, 1}}, {{25, 2, 3}, {45, 0, 1}}, {{3, 2, 3}, {45, 0, 1}}, {{0, 2, 3}, {45, 0, 1}}, {{25, 1, 3}, {45, 0, 1}}, {{25, 1, 2}, {45, 0, 1}}, {{3, 1, 3}, {45, 0, 1}}, {{3, 1, 2}, {45, 0, 1}}, {{0, 1, 3}, {45, 0, 1}}, {{0, 1, 2}, {45, 0, 1}}, {{25, 2, 3}, {36, 0, 1}}, {{3, 2, 3}, {36, 0, 1}}, {{0, 2, 3}, {36, 0, 1}}, {{25, 1, 3}, {36, 0, 1}}, {{25, 1, 2}, {36, 0, 1}}, {{3, 1, 3}, {36, 0, 1}}, {{3, 1, 2}, {36, 0, 1}}, {{0, 1, 3}, {36, 0, 1}}, {{0, 1, 2}, {36, 0, 1}}, {{25, 2, 3}, {23, 0, 1}}, {{3, 2, 3}, {23, 0, 1}}, {{0, 2, 3}, {23, 0, 1}}, {{25, 1, 3}, {23, 0, 1}}, {{25, 1, 2}, {23, 0, 1}}, {{3, 1, 3}, {23, 0, 1}}, {{3, 1, 2}, {23, 0, 1}}, {{0, 1, 3}, {23, 0, 1}}, {{0, 1, 2}, {23, 0, 1}}};

	coor final_coeff;

	// STEP 2 : computes the coefficient corresponding to the target monomial
	// For all trails, combine two products of size 2 to obtain a product of size 4
#pragma omp parallel for default(none) shared(target, list_trails, same_col_products, final_coeff, std::cout)
	for(const auto &[t0, t1]: list_trails) {
		if(!same_col_products[t0].empty() && !same_col_products[t1].empty()) {
			const coor cur_trail_product = multiply_maps_S6(same_col_products[t0], same_col_products[t1], target);
#pragma omp critical
			{
				final_coeff = add_coor(final_coeff, cur_trail_product);
			}
			cout << "|" << flush;
		}
	}
	cout << convert_monom_to_txt(final_coeff) << endl;
	return convert_monom_to_txt(final_coeff);
}