/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : rounds_1_to_4.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the functions needed to compute THE FIRST 4 ROUNDS of ASCON
 *           The computation is not exhaustive, it only computes the part of the
 *           ANF which is needed.
*/

#include "rounds_1_to_4.hpp"

using namespace std;

/*
 * For each row of the state s, prints the average length of a coordinate
 */
void print_len(const state &s, const uint &col, const string &name) {
	cout << name << " ";
	for(uint i = 0; i < 5; i++) {
		cout << s[i * 64 + col].size() << " ";
	}
	cout << endl;
}


/*
 *  Returns the addition of two coordinates.
 */
const coor add_coor(const coor &c1, const coor &c2) {
	coor c;
	set_symmetric_difference(c1.begin(), c1.end(), c2.begin(), c2.end(), std::inserter(c, c.begin()));
	return c;
}


/*
 * Returns the product of a monomial/monomial multiplication
 */
const monom mult_monom(const monom &m1, const monom &m2) {
	monom m = {0, 0, 0, 0, 0};
	for(uint i = 0; i < 5; i++)
		m[i] = m1[i] | m2[i];
	return m;
}


/*
 * Returns the product of a coordinate/coordinate multiplication.
 * The parameter condition_mult is a function  f: monomial -> Boolean
 * It is used to filter the resulting product: once the product of two monomials
 * is computed, we check if it is interesting or not for the next steps
 * (i.e. f(m1*m2) = true/false). All products of present monomials are computed
 *  during the multiplication of coordinates BUT only the interesting ones are
 *  stored in the resulting product.
 */
const coor mult_coor(const coor &c1, const coor &c2, const function<bool(const monom &)> &condition_mult) {
	coor c;
	for(const auto &x: c1) {
		for(const auto &y: c2) {
			const monom m = mult_monom(x, y);
			if(condition_mult(m)) {
				if(c.contains(m))
					c.erase(m);
				else
					c.insert(m);
			}
		}
	}
	return c;
}


/*
 * Returns true if the degree in public variable of monomial m is in the set of
 * degree degs.
 * This function is used to filter some multiplications of coordinates.
 */
bool cond_degree(const monom &m, const set<uint> &degs) {
	return degs.contains((uint) __builtin_popcountll(m[0]));
}


/*
 * ASCON Sbox function.
 * Input coordinates : x0 to x4
 * Output coordinates: y0 to y4
 * The Boolean parameter "quadrqtic" is used to indicate whether we need to
 * compute the whole Sbox layer or only the quadratic terms of the Sbox.
 * condition_mult is used to filter the resulting multiplications of coordinates.
 */
void sbox(const coor &x0, const coor &x1, const coor &x2, const coor &x3,
          const coor &x4, coor &y0, coor &y1, coor &y2, coor &y3, coor &y4,
          const bool &quadratic, const function<bool(const monom &)> &condition_mult) {
	const coor x2x1 = mult_coor(x2, x1, condition_mult);
	y2 = mult_coor(x4, x3, condition_mult);
	y3 = mult_coor(x0, add_coor(x3, x4), condition_mult);
	y4 = mult_coor(x1, add_coor(x4, x0), condition_mult);
	y1 = add_coor(mult_coor(add_coor(x2, x1), x3, condition_mult), x2x1);
	y0 = add_coor(x2x1, y4);

	if(!quadratic) {
		const coor x0_x1_x2_x3 = add_coor(add_coor(add_coor(x0, x1), x2), x3);
		const coor const_one = {{0,0,0,0,0}};
		y0 = add_coor(y0, x0_x1_x2_x3);
		y1 = add_coor(add_coor(y1, x0_x1_x2_x3), x4);
		y2 = add_coor(add_coor(add_coor(y2, x1), add_coor(x2, const_one)), x4);
		y3 = add_coor(add_coor(y3, x0_x1_x2_x3), x4);
		y4 = add_coor(add_coor(add_coor(y4, x1), x3), x4);

	}
}


/*
 * Returns the new state after applying the Sbox to the 64 columns of the state.
 */
const state sbox_state(const state &s, const bool &quadratic,
                       const function<bool(const monom &)> &condition_mult) {
	state new_state;

#pragma omp parallel for default(none) shared(s, new_state, std::cout, quadratic, condition_mult)
	for(uint i = 0; i < 64; i++) {
		sbox(s[i], s[i + 64], s[i + 128], s[i + 192], s[i + 256],
		     new_state[i], new_state[i + 64], new_state[i + 128],
		     new_state[i + 192], new_state[i + 256], quadratic, condition_mult);
	}
	return new_state;
}


/*
 * ASCON linear layer
 */
const state lin_layer(const state &s) {
	state new_state;
	const uint shifts[10] = {45, 36, 3, 25, 63, 58, 54, 47, 57, 23};

#pragma omp parallel for default(none) shared(s, new_state, shifts, std::cout)
	for(int j = 0; j < 64; ++j) {
		for(int i = 0; i < 5; ++i) {
			const uint cur = (i * 64) + j;
			new_state[cur] = add_coor(add_coor(s[cur], s[(i * 64) + ((j + shifts[i * 2]) % 64)]),
			                          s[(i * 64) + ((j + shifts[(i * 2) + 1]) % 64)]);

		}
		cout << "|" << flush;
	}
	cout << endl;
	return new_state;
}

/*
 * Converts a coordinate seen as a set of monomials into a polynomial whose
 * variables are v_i and coefficients are polynomials in 1, and a_i.
 * This corresponds to a usual F[x,y] = F[x][y] isomorphism.
 */
const poly_map convert_coor_to_poly_map(const coor &c) {
	map <uint64_t, coor> m;

	for(const auto &x: c) {
		uint64_t cur_monom = x[0];
		if(m[cur_monom].contains(x))
			m[cur_monom].erase(x);
		else
			m[cur_monom].insert(x);
	}
	return m;
}


/*
 * Returns the state after the fourth linear layer from a given initial state.
 */
const state build_state_l4(const state &start) {

	// Filter function : does not filter anything
	const function<bool(const monom &)> f_s1 = [](const monom &m) { return true; };
	const state s1 = sbox_state(start, false, f_s1); // true Sbox
	print_len(s1, 0, "s1");
	const state l1 = lin_layer(s1);
	print_len(l1, 0, "l1");

	const set<uint> deg_s2 = {2};
	// Filter function : only product of degree 2
	const function<bool(const monom &)> f_s2 = [deg_s2](const monom &m) { return cond_degree(m, deg_s2); };
	const state s2 = sbox_state(l1, true, f_s2); // BEWARE, highest-degree terms so quad is true
	print_len(s2, 0, "s2");
	const state l2 = lin_layer(s2);
	print_len(l2, 0, "l2");

	const set<uint> deg_s3 = {4};
	// Filter function : only product of degree 4
	const function<bool(const monom &)> f_s3 = [deg_s3](const monom &m) { return cond_degree(m, deg_s3); };
	const state s3 = sbox_state(l2, true, f_s3); // quadratic Sbox
	print_len(s3, 0, "s3");
	const state l3 = lin_layer(s3);
	print_len(l3, 0, "l3");

	const set<uint> deg_s4 = {8};
	// Filter function : only product of degree 8
	const function<bool(const monom &)> f_s4 = [deg_s4](const monom &m) { return cond_degree(m, deg_s4); };
	const state s4 = sbox_state(l3, true, f_s4); // quadratic Sbox
	print_len(s4, 0, "s4");
	const state l4 = lin_layer(s4);
	print_len(l4, 0, "l4");

	return l4;
}


/*
 * Converts the state into an array of poly_maps
 */
const array<poly_map, 320> convert_l4(const state &l4) {
	cout << "conversion..." << endl;

	array<poly_map, 320> l4_converted;
#pragma omp parallel for default(none) shared(l4_converted, l4, std::cout)
	for(uint i = 0; i < 320; i++) {
		l4_converted[i] = convert_coor_to_poly_map(l4[i]);
		cout << "|" << flush;
	}
	cout << endl;
	return l4_converted;
}


/*
 * Returns the state after the fourth linear layer as an array of poly_map
 * from a given initial state.
 */
const array<poly_map, 320> get_l4(const state &start) {
	return convert_l4(build_state_l4(start));
}
