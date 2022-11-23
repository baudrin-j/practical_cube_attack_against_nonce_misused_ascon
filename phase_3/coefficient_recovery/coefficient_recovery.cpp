/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : coefficient_recovery.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Contains the main function for the recovery of coefficients of
 * degree-31 monomials after the sixth S-box layer under the assumption that
 * vectors a and e are already fully-recovered.
*/

#include "coefficient_recovery.hpp"

using namespace std;
using namespace chrono;

// Random 64-bit word generator
uint64_t random_monom() {
	uint64_t m = 0;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> distrib(0, 1);
	for(uint i = 0; i < 64; i++) {
		uint t = distrib(gen);
		if(t)
			m |= ((uint64_t) 1) << i;
	}
	return m;
}

/*
 * Initialize a state with cube variables contained in "cube" only, all the
 * a_i and e_i VALUES, and all b_i and c_i VARIABLES.
 */
const state initialize_state(set<uint> cube, set<uint> list_a, set<uint> list_e_0) {
	state start;
	for(uint i = 0; i < 5; i++) {
		for(uint j = 0; j < 64; j++) {
			// ROW 0 -- v_i inserted only for i in cube
			if((i == 0 && cube.count(j))) {
				monom v = {0, 0, 0, 0, 0};
				v[0] = ((uint64_t) 1) << (63 - j);
				start[j].insert(v);
			}
			// ROW 1 -- a, a_i inserted as constants
			if((i == 1 && list_a.count(j))) {
				monom one = {0, 0, 0, 0, 0};
				start[64 + j].insert(one);
			}
			// ROW 2, 3 -- b/c, b_i and c_i inserted as variables
			if(i == 2 || i == 3) {
				monom bc = {0, 0, 0, 0, 0};
				bc[i] = ((uint64_t) 1) << (63 - j);
				start[i * 64 + j].insert(bc);
			}
			// ROW 4 -- d = c + (e + 1), c_i inserted as variables, e_i as constant
			if(i == 4) {
				monom c = {0, 0, 0, 0, 0};
				c[3] = ((uint64_t) 1) << (63 - j);
				start[i * 64 + j].insert(c);
				if(list_e_0.count(j)) {
					monom one = {0, 0, 0, 0, 0};
					start[i * 64 + j].insert(one);
				}
			}
		}
	}
	return start;
}


int main() {
	omp_set_num_threads(8);

	// Random a, e with uniformly distributed a_i, e_i bits
	set<uint> list_a; // List of i such that a_i = 1
	set<uint> list_e_1; // List of i such that e_i = 1
	set<uint> list_e_0; // List of i such that e_i = 0

	uint64_t a = 0; // Mask corresponding to a
	uint64_t e = 0; // Mask corresponding to e

	for(uint i = 0; i < 64; i++) {
		if(random_monom() % 2) {
			list_a.insert(i);
			a |= (((uint64_t) 1) << (63 - i));
		}
		if(random_monom() % 2) {
			list_e_1.insert(i);
			e |= (((uint64_t) 1) << (63 - i));
		} else
			list_e_0.insert(i);
	}

	/*
	 * Build random cubes made as follows:
	 * Add ``as much v_i as possible such that e_i = 0'' while assuring a minimal
	 * number of v_i such that e_i = 1
	 */
	const uint nb_cubes = 3;
	array <set<uint>, nb_cubes> cubes; // cubes as lists of integers
	array <uint64_t, nb_cubes> targets; // cubes as binary masks
	const uint nb_zeros = 28; // maximum number of v_i such that e_i = 0. nb_zeros < 30
	for(uint i = 0; i < nb_cubes; i++) {
		// As much v_i as possible such that e_i = 0 while respecting the max nb.
		if(list_e_0.size() > nb_zeros) {
			while(cubes[i].size() != nb_zeros) {
				cubes[i].insert(*next(list_e_0.begin(), random_monom() % list_e_0.size()));
			}
		}
		else
			cubes[i].insert(list_e_0.begin(), list_e_0.end());

		// Add at least one var v_i such that e_i = 1, more if necessary
		while(cubes[i].size() != 31)
			cubes[i].insert(*next(list_e_1.begin(), random_monom() % list_e_1.size()));

		// Fill the mask
		for(auto &x: cubes[i])
			targets[i] |= ((uint64_t) 1) << (63 - x);
	}

	// Quick verification, we need distinct cubes
	for(uint i = 0; i < nb_cubes; i++) {
		for(uint j = i + 1; j < nb_cubes; j++) {
			if(cubes[i] == cubes[j]) {
				cout << "Bad choice of cubes" << endl;
				return 1;
			}
		}
	}

	/*
	 * Save the parameters in a file
	 * Format:
	 *          - row 0 : a given in hex
	 *          - row 1 : e given in hex
	 *          - rows 2 to 2 + nb_cubes: cubes given in hex
	 */
	ofstream parameters("../results/parameters.txt");
	parameters << std::hex << a << endl << e << endl;
	for(uint i = 0; i < nb_cubes; i++) {
		parameters << targets[i] << endl;
	}
	parameters.close();

	for(uint k = 0; k < nb_cubes; k++) {
		/*
		 * STEP 1: For each cube, initialize a state with:
		 *          - the necessary variables (all b_i, all c_i, and 31 v_i)
		 *          - the values of a_i and e_i
		*/
		const state start = initialize_state(cubes[k], list_a, list_e_0);

		// STEP 2: Compute all the terms of deg 7 or 8 after L4
		const array<poly_map, 320 > l4 = get_l4(start);

		// STEP 3: Compute the coefficients of the targeted cube of degree 31 after S6.
		auto start_step3 = high_resolution_clock::now();

		// FOR LIMITED MEMORY USAGE : compute the coefficients one by one
		for(int i = 0; i < 64; ++i) {
			const auto start_col = high_resolution_clock::now();
			cout << "Col" + to_string(i) + "..." << endl;
			ofstream f("../results/polynomials_cube_" + to_string(k) + ".txt", fstream::out | fstream::app);
			f << coefficient_recovery(i, l4, targets[k]);
			f.close();
			const auto stop_col = high_resolution_clock::now();
			const auto duration_col = duration_cast<seconds>(stop_col - start_col);
			cout << k << ", col" << i << " done in " + to_string(duration_col.count()) + "secs" << endl;
		}

		// FOR ENVIRONMENT NOT LIMITED BY MEMORY : compute all the coefficients in parallel
		// coefficient_recovery_all_polys(l4, targets[k], "../results/polynomials_cube_" + to_string(k) + ".txt");

		auto stop_step3 = high_resolution_clock::now();
		auto duration_step3 = duration_cast<seconds>(stop_step3 - start_step3);
		cout << "\n\n S5-L5-S6 in " + to_string(duration_step3.count()) + "secs\n";
	}
	return 0;
}
