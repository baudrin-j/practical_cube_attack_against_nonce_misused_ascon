/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : coefficient_recovery.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Contains the main function for the recovery of coefficients of
 * degree-32 monomials after the sixth S-box layer under the assumption that all
 * the bits of e and about half the bits of a are already recovered.
*/

#include "coefficient_recovery.hpp"

using namespace std;
using namespace chrono;


//// Random 64-bit word generator
//uint64_t random_monom() {
//	uint64_t m = 0;
//	std::random_device rd;  //Will be used to obtain a seed for the random number engine
//	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//	std::uniform_int_distribution<> distrib(0, 1);
//	for(uint i = 0; i < 64; i++) {
//		uint t = distrib(gen);
//		if(t)
//			m |= ((uint64_t) 1) << i;
//	}
//	return m;
//}


set<uint> select_cube(uint *nb_unknowns, uint64_t *target, uint64_t a, uint64_t e, set<uint> &list_e_0, set<uint> &list_e_1, set<uint> &list_a_recovered){
	set<uint> cube;
	*target = (uint64_t) 0; // Maks corresponding to cubes
	const uint nb_zeros = 29; // CAN BE MODIFIED

	// Choice of cube : at least one var v_i such that e_i = 1, as many v_i as possible such that e_i = 0
	if(list_e_0.size() > nb_zeros) { // add a random subset of size nb_zeros
		while(cube.size() != nb_zeros) {
			cube.insert(*next(list_e_0.begin(), random_monom() % list_e_0.size()));
		}
		cout << "NB_ZEROS == " << nb_zeros <<  endl;
	}
	else { // add as much variables st e_i = 0 as possible
		cube.insert(list_e_0.begin(), list_e_0.end());
		cout << "NB_ZEROS == " << list_e_0.size() << endl;
	}

	while(cube.size() != 32) {// Add at least one var v_i such that e_i = 1, more if necessary
		uint r = random_monom() % list_e_1.size();
		uint tmp_index = *next(list_e_1.begin(), r);
		if(!list_a_recovered.count(tmp_index)) {
			cube.insert(tmp_index);
			(*nb_unknowns)++;
		}
	}
	cout << "NB unknowns: " << *nb_unknowns << endl;

	// Fill the mask
	for(auto &x: cube)
		*target |= ((uint64_t) 1) << (63 - x);

/*
 * Save the parameters in a file
 * Format:
 *          - row 0 : a given in hex
 *          - row 1 : e given in hex
 *          - rows 2 : current cube in hex
 */
	ofstream parameters("results/parameters.txt");
	parameters << std::hex << a << endl << e << endl;
	parameters << *target << endl;
	parameters.close();

	return cube;
}

/*
 * Initialize a state with cube variables contained in "cube" only, as well as
 * all the a_i and e_i with i in "cube" only.
 * Depending on their status, a_i and e_i are input as constant (if they have
 * already been recovered during phase 1, or since the beginning of this phase),
 * or as variable.
 */
state initialize_state(const set<uint> &cube, const set<uint> &list_a,
					   const set<uint> &list_e_0, const set<uint> &list_a_recovered_1){
	state start;

	for(uint i = 0; i < 5; i++) {
		for(uint j = 0; j < 64; j++) {
			//ROW 0 -- v
			if((i == 0 && cube.count(j))) {
				monom v = {0, 0, 0, 0, 0};
				v[0] = ((uint64_t) 1) << (63 - j);
				start[j].insert(v);
			}
			//ROW 1 -- a
			if((i == 1 && cube.count(j))){
				if((list_e_0.count(j) && list_a.count(j)) || list_a_recovered_1.count(j)) {
					// if e = 0 inserted as value
					monom one = {0, 0, 0, 0, 0};
					start[64 + j].insert(one);
				}
				else if(!list_e_0.count(j)) { // if e = 1, inserted as variable
					monom a = {0, 0, 0, 0, 0};
					a[1] = ((uint64_t) 1) << (63 - j);
					start[64 + j].insert(a);
				}
			}
			//ROW 3,4 -- c inserted as variables
			if(i == 3 || i == 4) {
				monom c = {0, 0, 0, 0, 0};
				c[3] = ((uint64_t) 1) << (63 - j);
				start[i * 64 + j].insert(c);
				if(i == 4 && list_e_0.count(j)) {
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
	uint max_tries = 15;

	//STEP 0 : Initialization of capacity rows a & e
	set<uint> list_a; // List of i such that a_i = 1
	set<uint> list_e_1; // List of i such that e_i = 1
	set<uint> list_e_0; // List of i such that e_i = 0
	set<uint> list_a_recovered; // List of recovered a_i on the fly
	set<uint> list_a_recovered_1; // List of recovered a_i on the fly such that a_i = 1

	uint64_t a = 0; // Mask corresponding to a
	uint64_t e = 0; // Mask corresponding to e

	// Random a, e with uniformly distributed a_i, e_i bits
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
	cout << "Total number of unknowns a_i:" << list_e_1.size() << endl;

	// Loop until there is no more bit to recover or the max number of tries is reached.
	uint tries = 0;
	while(list_a_recovered != list_e_1 && tries < max_tries) {
		tries++;
		uint nb_unknowns = 0;
		uint64_t target = 0;
		set<uint> cube = select_cube(&nb_unknowns, &target, a, e, list_e_0, list_e_1, list_a_recovered);

		/*
		* STEP 1: For each cube, initialize a state with:
        *          - the necessary variables (all b_i, all c_i, and 31 v_i)
		*          - the values of a_i and e_i
		*/
		state start = initialize_state(cube, list_a, list_e_0, list_a_recovered_1);

		// STEP 2: Compute all the terms of deg 8 after L4
		const array<poly_map , 320> l4 = get_l4(start);

		// STEP 3: Compute the coefficients of the targeted cube of degree 32 after S6.
		uint count_non_constant = 0;
		for(int i = 0; i < 64; i++) {
			string s = coefficient_recovery(i, l4, target);

			ofstream f;
			if(i)
				f.open("results/polynomials.txt",fstream::out | fstream::app);
			else
				f.open("results/polynomials.txt");
			f << s << endl;
			f.close();

			if(s != "0" && s != "1")
				count_non_constant++;

			// Automatic stop if the number of non-constant equations is twice
			// the number of unknowns
			if(count_non_constant > (2 * nb_unknowns))
				break;
		}

		// STEP 4 : Compute the corresponding cube-sum
		cout << "values recovery..." << endl;
		cube_sum_given_cubes_given_a_e("results/parameters.txt", "results/cube_sum_vectors.txt");

		// STEP 5 : From the polynomials and the values, build the system and
		// solve it.
		cout << "system solving..." << endl;
		system("zsh script.run"); // SHOULD BE MODIFIED IF ANOTHER SHELL IS USED
		cout << nb_unknowns << "|| " << list_a_recovered.size() << endl;

		// STEP 6 : Take the newly recovered values into account
		ifstream recovered_a("results/recovered_a.txt");
		string line;
		while(getline (recovered_a,line)) {
			string var = line.substr(1, line.find(" = ") - 1);
			list_a_recovered.insert(stoi(var));
			uint value = stoi(line.substr(var.size() + 4, line.size() - (var.size() + 4)));
			if(value)
				list_a_recovered_1.insert(stoi(var));
		}
		recovered_a.close();
		cout << nb_unknowns << "|| " << list_a_recovered.size() << endl;
	}
	return 0;
}