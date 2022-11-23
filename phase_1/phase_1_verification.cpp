/*
 * Filename : phase_1_verification.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Main file for the verification of the first phase.
*/
#include "cube_computation.h"

using namespace std;
using namespace std::chrono;


/*
 * Launches some trials for any of the two cubes introduced in our paper.
 * The results are stored in a folder called "results" (PLEASE CREATE THE FOLDER BEFORE)
 * The program takes as input two parameters:
 * - a header for the result files which will be named
 *   {header}_cube_a_{0,1}_e_{0,1}.txt. It can for instance be a subfolder name.
 * - The second parameter is the index corresponding to the cube:
 *      - 0 stands for cube x^v from our paper.
 *      - any other integer stands for cube x^w from our paper.
 * - As the last linear layer can be inverted, it is omitted.
 * - Variable "nb_tries" can be increased if a higher number of trials are needed.
 */
int main(int argc, char *argv[]){
	if(argc != 3)
		return 1;

	omp_set_num_threads(8);
	uint rounds = 6;
	bool last_lin = false;
	bool cst = true;
	uint nb_tries = 10; // can be modified according to the needs
	const string header = argv[1];
	const uint cube_index = stoi(argv[2]);

	vector<uint> cube;
	if(cube_index == 0)
		cube = {0, 1, 4, 5, 6, 8, 14, 15, 16, 26, 27, 30, 34, 37, 38, 48, 49, 50, 56, 58, 59, 60, 63, 17, 35, 40, 46, 55, 9, 12, 18, 19};
	else
		cube = {0, 1, 4, 5, 6, 8, 14, 15, 16, 26, 27, 30, 34, 37, 38, 48, 49, 50, 56, 58, 59, 60, 63, 17, 35, 40, 46, 55, 7, 24, 41, 43};

	/*
	 * For "nb_tries" random capacities, the cube-sum corresponding to x^v or x^w
	 * is computed. Then the cube-sum vector is stored as a hex string in one
	 * of the four output files (a_0_e_0, a_0_e_1, a_1_e_0, a_1_e_1), depending
	 * on the values of the most significant bits of a and e.
	 * From our experimentations, only a single file contains all-zero output
	 * vectors.
	 * */

	for(uint i = 0; i < nb_tries; i++) {
		auto start = high_resolution_clock::now();

		// New random inner state
		uint64_t state[5] = {0,0,0,0,0};
		for(uint j = 1; j < 5; j++)
			state[j] = random_monom();

		uint e = ((~(state[3] ^ state[4])) >> 63) & 1;
		uint a = (state[1] >> 63) & 1;

		cube_sum(state, rounds, cube, last_lin, cst);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);

		// File saving
		ofstream f("results/" + header + "_cube_" + to_string(cube_index) + "_a_" + to_string(a) + "_e_" + to_string(e) + ".txt", fstream::out | fstream::app);
		f << std::hex << state[0] << endl;
		f.close();

		// Quick overview of the current result
		cout << i << " Time: " << duration.count() << " | a: " << a <<  " | e: " << e << " | w:" <<  __builtin_popcountll(state[0]) << endl;
	}
	return 0;
}
