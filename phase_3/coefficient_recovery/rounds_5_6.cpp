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
string convert_monom_to_txt(const coefficient &m) {
	string s;
	bool first = true;
	if(m[0]) { // Handles the constant, row 0 of the coefficient
		s = "1";
		first = false;
	}

	for(uint i = 1; i < 4; i++) { // Handles rows 1, 2 and 3
		if(m[i]) {
			for(uint j = 0; j < 64; j++) {
				if((m[i] >> (63 - j)) & 1) {
					if(first)
						first = false;
					else
						s += " + ";

					if(i == 1)
						s += "b" + to_string(j) + "*c" + to_string(j);
					else if(i == 2)
						s += "b" + to_string(j);
					else
						s += "c" + to_string(j);
				}
			}
		}
	}

	if(first) // Handles the case of the null coefficient
		s = "0\n";
	else
		s += "\n";

	return s;
}

/*
 * Computes the addition of two coefficients represented as an array of 4 uint64_t.
 * Corresponds to += operator for coefficients.
 */
void add_coeff(coefficient &coeff, const coefficient &coeff_to_add) {
	for(uint i = 0; i < 4; i++) {
		coeff[i] = (coeff[i] ^ coeff_to_add[i]);
	}
}


/*
 * Computes a partial multiplication between two poly_map and returns a poly_map.
 * It is expected that c1 and c2 are two poly_maps with terms of degree 7 or 8, as output by get_l4.
 * It only returns the terms of degree 15 or more that appears in the product.
 *
 * This function corresponds to the computation of the interesting terms during S5.
 */
const poly_map multiply_maps_S5(const poly_map &c1, const poly_map &c2) {
	poly_map prod; // Output product

	// Double for loop to compute the product, restricted by the condition deg(tmp_monom) >= 15
	for(const auto &[monom1, coeff1]: c1) {
		const bool add_coeff1 = (((uint) __builtin_popcountll(monom1)) == 7);

		for(const auto &[monom2, coeff2]: c2) {
			const uint64_t tmp_monom = (monom1 | monom2); // Multiplication of two monomials is an OR
			if(((uint) __builtin_popcountll(tmp_monom)) >= 15) {
				// Initial value in the map if tmp_monom was not present yet
				if(!prod.contains(tmp_monom))
					prod[tmp_monom] = {(uint64_t) 0, (uint64_t) 0, (uint64_t) 0, (uint64_t) 0};

				// Addition of coeff1*coeff2 to the coefficient of tmp_monom.
				// In our case coeff1 and/or coeff2 is equal to 1 (depending on the degree of monom1/2)
				// So the multiplication vanished, but we need to select the (possibly) non-constant coefficient.
				if(add_coeff1)
					add_coeff(prod[tmp_monom], coeff1);
				else
					add_coeff(prod[tmp_monom], coeff2);
			}
		}
	}
	return prod;
}


/*
 * Computes a partial multiplication between two poly_map and returns a coefficient.
 * It is expected that c1 and c2 are two poly_maps with terms of degree 15 or 16, as output by multiply_maps_S5.
 * It only returns the coefficient (that appears in the product c1*c2) corresponding to the target monomial given as input.
 * It is expected that target is of degree (hamming weight) 31.
 *
 * This function corresponds to the computation of a coefficient of a monomial of degree 31 after S6.
 */
coefficient multiply_maps_S6(const poly_map &c1, const poly_map &c2, const uint64_t &target) {
	coefficient prod = {0, 0, 0, 0}; // Output coefficient

	// Select the smallest list to be browsed
	const poly_map * first = &c1;
	const poly_map * second = &c2;
	if(c2.size() < c1.size()) {
		first = &c2;
		second = &c1;
	}

	for(const auto &[monom1, coeff1]: (*first)) { // Loop over the smallest list
		if(coeff1[0] || coeff1[1] || coeff1[2] || coeff1[3]) { // If monom1 actually appears
			const bool monom1_subleading = (((uint) __builtin_popcountll(monom1)) == 15);
			const uint64_t complement = ((~monom1) & target);

			if((*second).contains(complement)) { // Look for the complementary monom in the second list
				const coefficient &coeff2 = (*second).at(complement);

				if(coeff2[0] || coeff2[1] || coeff2[2] || coeff2[3]) { // If complement actually appears
					if(monom1_subleading)
						add_coeff(prod, coeff1);
					else
						add_coeff(prod, coeff2);
				}
			}

			if(!monom1_subleading) { // If deg(monom1) == 16, look for the possible covering of the target by two monoms of deg 16
				for(uint i = 0; i < 64; i++) {
					if((monom1 >> i) & 1) {
						const uint64_t covering = (complement | (((uint64_t) 1) << i));
						if((*second).contains(covering) && (*second).at(covering)[0]) // If the covering monomial is present and its coeff is equal to 1 (= non-null)
								add_coeff(prod, {(uint64_t) 1, (uint64_t) 0, (uint64_t) 0, (uint64_t) 0});
					}
				}
			}
		}
	}
	return prod;
}


/*
 * Computes the coefficient of a monomial of degree 31 in a single coordinate c_{0,x} after S6.
 * The coefficient is output as a string.
 *
 * - col is the index of the coordinate in which we are looking for. (0 <= col <= 63)
 * - target is the monomial we are targeting. (target is a word of size 64 and hamming weight 31)
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

	coefficient final_coeff = {(uint64_t) 0, (uint64_t) 0, (uint64_t) 0, (uint64_t) 0};

	// STEP 2 : computes the coefficient corresponding to the target monomial
	// For all trails, combine two products of size 2 to obtain a product of size 4
#pragma omp parallel for default(none) shared(target, list_trails, same_col_products, final_coeff, std::cout)
	for(const auto &[t0, t1]: list_trails) {
		if(!same_col_products[t0].empty() && !same_col_products[t1].empty()) {
			const coefficient cur_trail_product = multiply_maps_S6(same_col_products[t0], same_col_products[t1], target);
			add_coeff(final_coeff, cur_trail_product);
			cout << "|" << flush;
		}
	}
	cout << endl << "final length of the polynomial :" << ((uint) __builtin_popcountll(final_coeff[0]) + (uint) __builtin_popcountll(final_coeff[1]) + (uint) __builtin_popcountll(final_coeff[2]) + (uint) __builtin_popcountll(final_coeff[3])) << endl;
	return convert_monom_to_txt(final_coeff);
}


/*
 * Computes the 64 coefficients of a monomial of degree 31 present on row 0 after S6.
 * The coefficients are stored in a text file.
 *
 * - target is the monomial we are targeting. (target is a word of size 64 and hamming weight 31)
 * - l4 is the state after l4, initialized with only the necessary variables.
 *  It is expected that l4 has been initialized through get_l4 first.
 * - filename contains the outputfile location.
 */
void coefficient_recovery_all_polys(const array<poly_map, 320> &l4, const uint64_t &target, const string &filename) {
	const auto start_s5 = high_resolution_clock::now();
	cout << "S5/L5..." << endl;

	// List of the necessary products occurring during S5 in each column.
	// Each generic_size_2_products corresponds to the indexes of two rows multiplied through ASCON S-box.
	// The list is not exhaustive as not all products appear in the 1.5-round trails leading to the first row.
	const vector <generic_size_2_products> list_products = {{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {3, 4}};

	// Table mapping a product to its actual polynomial
	array<array<poly_map, 6>, 64> same_col_products;

	// Step 1 : for each col and for each product, compute the product and store it in the table.
#pragma omp parallel for default(none) shared(l4, list_products, same_col_products, std::cout)
	for(uint j = 0; j < 64; j++) {
		const auto start_col = high_resolution_clock::now();
		for(uint i = 0; i < list_products.size(); i++) {
			const auto &[y1, y2] = list_products[i];
			same_col_products[j][i] = multiply_maps_S5(l4[y1 * 64 + j], l4[y2 * 64 + j]);
		}

		const auto stop_col = high_resolution_clock::now();
		const auto duration_col = duration_cast<seconds>(stop_col - start_col);
		cout << "Col " + to_string(j) + " - done in " + to_string(duration_col.count()) + "secs" << endl;
	}

	const auto stop_s5 = high_resolution_clock::now();
	const auto duration_s5 = duration_cast<seconds>(stop_s5 - start_s5);
	cout << "S5-L5 done in " + to_string(duration_s5.count()) + "secs.\nS6...";

	// List of trails of size 4 leading to coordinate c_{0,0} through 1.5 round
	const vector<trails> list_trails = {{{25, 2, 3}, {63, 3, 4}}, {{25, 2, 3}, {58, 3, 4}}, {{0, 3, 4}, {25, 2, 3}}, {{3, 2, 3}, {63, 3, 4}}, {{58, 3, 4}, {3, 2, 3}}, {{0, 3, 4}, {3, 2, 3}}, {{0, 2, 3}, {63, 3, 4}}, {{0, 2, 3}, {58, 3, 4}}, {{57, 1, 4}, {25, 2, 3}}, {{57, 1, 4}, {3, 2, 3}}, {{0, 2, 3}, {57, 1, 4}}, {{25, 2, 3}, {45, 1, 4}}, {{25, 2, 3}, {45, 1, 2}}, {{3, 2, 3}, {45, 1, 4}}, {{3, 2, 3}, {45, 1, 2}}, {{0, 2, 3}, {45, 1, 4}}, {{0, 2, 3}, {45, 1, 2}}, {{25, 2, 3}, {36, 1, 4}}, {{25, 2, 3}, {36, 1, 2}}, {{3, 2, 3}, {36, 1, 4}}, {{3, 2, 3}, {36, 1, 2}}, {{0, 2, 3}, {36, 1, 4}}, {{0, 2, 3}, {36, 1, 2}}, {{25, 1, 3}, {63, 3, 4}}, {{25, 1, 3}, {58, 3, 4}}, {{0, 3, 4}, {25, 1, 3}}, {{25, 1, 2}, {63, 3, 4}}, {{25, 1, 2}, {58, 3, 4}}, {{0, 3, 4}, {25, 1, 2}}, {{25, 1, 3}, {57, 1, 4}}, {{25, 1, 2}, {57, 1, 4}}, {{25, 1, 3}, {45, 1, 4}}, {{25, 1, 3}, {45, 1, 2}}, {{25, 1, 2}, {45, 1, 4}}, {{25, 1, 2}, {45, 1, 2}}, {{25, 1, 3}, {36, 1, 4}}, {{25, 1, 3}, {36, 1, 2}}, {{25, 1, 2}, {36, 1, 4}}, {{25, 1, 2}, {36, 1, 2}}, {{25, 2, 3}, {23, 1, 4}}, {{3, 2, 3}, {23, 1, 4}}, {{0, 2, 3}, {23, 1, 4}}, {{25, 1, 3}, {23, 1, 4}}, {{25, 1, 2}, {23, 1, 4}}, {{3, 1, 3}, {63, 3, 4}}, {{58, 3, 4}, {3, 1, 3}}, {{0, 3, 4}, {3, 1, 3}}, {{3, 1, 2}, {63, 3, 4}}, {{58, 3, 4}, {3, 1, 2}}, {{0, 3, 4}, {3, 1, 2}}, {{57, 1, 4}, {3, 1, 3}}, {{57, 1, 4}, {3, 1, 2}}, {{3, 1, 3}, {45, 1, 4}}, {{3, 1, 3}, {45, 1, 2}}, {{3, 1, 2}, {45, 1, 4}}, {{3, 1, 2}, {45, 1, 2}}, {{3, 1, 3}, {36, 1, 4}}, {{3, 1, 3}, {36, 1, 2}}, {{3, 1, 2}, {36, 1, 4}}, {{3, 1, 2}, {36, 1, 2}}, {{3, 1, 3}, {23, 1, 4}}, {{3, 1, 2}, {23, 1, 4}}, {{0, 1, 3}, {63, 3, 4}}, {{0, 1, 3}, {58, 3, 4}}, {{0, 1, 2}, {63, 3, 4}}, {{0, 1, 2}, {58, 3, 4}}, {{0, 1, 2}, {0, 3, 4}}, {{0, 1, 2}, {25, 2, 3}}, {{0, 1, 2}, {3, 2, 3}}, {{0, 1, 3}, {57, 1, 4}}, {{0, 1, 2}, {57, 1, 4}}, {{0, 1, 3}, {45, 1, 4}}, {{0, 1, 3}, {45, 1, 2}}, {{0, 1, 2}, {45, 1, 4}}, {{0, 1, 2}, {45, 1, 2}}, {{0, 1, 3}, {36, 1, 4}}, {{0, 1, 3}, {36, 1, 2}}, {{0, 1, 2}, {36, 1, 4}}, {{0, 1, 2}, {36, 1, 2}}, {{0, 1, 2}, {25, 1, 3}}, {{0, 1, 2}, {25, 1, 2}}, {{0, 1, 3}, {23, 1, 4}}, {{0, 1, 2}, {23, 1, 4}}, {{0, 1, 2}, {3, 1, 3}}, {{0, 1, 2}, {3, 1, 2}}, {{57, 0, 1}, {25, 2, 3}}, {{57, 0, 1}, {3, 2, 3}}, {{0, 2, 3}, {57, 0, 1}}, {{57, 0, 1}, {25, 1, 3}}, {{57, 0, 1}, {25, 1, 2}}, {{57, 0, 1}, {3, 1, 3}}, {{57, 0, 1}, {3, 1, 2}}, {{0, 1, 3}, {57, 0, 1}}, {{0, 1, 2}, {57, 0, 1}}, {{25, 2, 3}, {45, 0, 1}}, {{3, 2, 3}, {45, 0, 1}}, {{0, 2, 3}, {45, 0, 1}}, {{25, 1, 3}, {45, 0, 1}}, {{25, 1, 2}, {45, 0, 1}}, {{3, 1, 3}, {45, 0, 1}}, {{3, 1, 2}, {45, 0, 1}}, {{0, 1, 3}, {45, 0, 1}}, {{0, 1, 2}, {45, 0, 1}}, {{25, 2, 3}, {36, 0, 1}}, {{3, 2, 3}, {36, 0, 1}}, {{0, 2, 3}, {36, 0, 1}}, {{25, 1, 3}, {36, 0, 1}}, {{25, 1, 2}, {36, 0, 1}}, {{3, 1, 3}, {36, 0, 1}}, {{3, 1, 2}, {36, 0, 1}}, {{0, 1, 3}, {36, 0, 1}}, {{0, 1, 2}, {36, 0, 1}}, {{25, 2, 3}, {23, 0, 1}}, {{3, 2, 3}, {23, 0, 1}}, {{0, 2, 3}, {23, 0, 1}}, {{25, 1, 3}, {23, 0, 1}}, {{25, 1, 2}, {23, 0, 1}}, {{3, 1, 3}, {23, 0, 1}}, {{3, 1, 2}, {23, 0, 1}}, {{0, 1, 3}, {23, 0, 1}}, {{0, 1, 2}, {23, 0, 1}}};

	// List of output coefficients
	array<coefficient, 64> final_coors;
	for(uint i = 0; i < 64; i++) {
		for(uint j = 0; j < 4; j++) {
			final_coors[i][j] = (uint64_t) 0;
		}
	}

	ofstream f(filename, fstream::out | fstream::app);

	// Step 2: for each col and each trails of size 4 leading to the selected coordinate,
	// computes the current product of size 4 (as the product of two products of size 2)
	// and sums it with the partial coefficient
	for(uint i = 0; i < 64; i++) {
		const auto start_col = high_resolution_clock::now();

#pragma omp parallel for default(none) shared(target, list_trails, list_products, same_col_products, final_coors, std::cout, i)
		for(const auto &[t0, t1]: list_trails) {
			const auto &[x1, y1, y2] = t0;
			const uint prod1 = distance(list_products.begin(), find(list_products.begin(), list_products.end(), make_tuple(y1, y2)));
			const auto &[x2, y3, y4] = t1;
			const uint prod2 = distance(list_products.begin(), find(list_products.begin(), list_products.end(), make_tuple(y3, y4)));
			const auto &c1 = same_col_products[(x1 + i) % 64][prod1];
			const auto &c2 = same_col_products[(x2 + i) % 64][prod2];
			if(!c1.empty() && !c2.empty()) {
				const coefficient tmp_coeff = multiply_maps_S6(c1, c2, target);
#pragma omp critical
				{
					add_coeff(final_coors[i], tmp_coeff);
				}
			}
		}
		const auto stop_col = high_resolution_clock::now();
		const auto duration_col = duration_cast<seconds>(stop_col - start_col);
		cout << endl << "Poly " + to_string(i) + "in " + to_string(duration_col.count()) + "secs" << endl;
		f << convert_monom_to_txt(final_coors[i]);
	}
	f.close();
}