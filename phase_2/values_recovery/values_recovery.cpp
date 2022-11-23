/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : values_recovery.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Recovery of the cube-sum vector corresponding to a "parameter" file
 *           which was output by coefficient_recovery
*/

#include "values_recovery.h"

using namespace std;
using namespace std::chrono;

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
 * Given an input "parameters" file containing the already-recovered vectors a
 * and e, as well as a cube of size 32 this functions :
 *  1- first selects pseudo-random values for vectors b and c;
 *  2- then, computes the cube-sum vector corresponding to the cube stored
 *  in the input file;
 *  3- and finally, outputs a file containing the cube-sum vector.
 */
void cube_sum_given_cubes_given_a_e(const string &inputfilename, const string &outputfilename){
	uint rounds = 6;
	bool last_lin = false;
	bool cst = false;

	ifstream inputfile(inputfilename);
	string line;
	vector<uint64_t> lines;
	while(getline (inputfile,line))
		lines.insert(lines.end(), strtoull(line.c_str(), NULL, 16));
	inputfile.close();

	uint64_t a = lines[0];
	uint64_t e = lines[1];
	uint64_t cube_int = lines[2];

	// Random values are used for b and c, because the cube-sum value is
	// independent of them.
	uint64_t b = random_monom();
	uint64_t c = random_monom();

	vector<uint> cube;
	for(uint j = 0; j < 64; j++) {
		if((cube_int >> (63 - j)) & 1)
			cube.insert(cube.end(), j);
	}

	uint64_t state[5] = {0,0,0,0,0};
	state[1] = a;
	state[2] = b;
	state[3] = c;
	state[4] = ~(c ^ e);

	cube_sum(state, rounds, cube, last_lin, cst);
	ofstream outputfile(outputfilename);
	outputfile << std::hex << state[0] << endl;
	outputfile.close();
}

/* int main(){

	omp_set_num_threads(8);
	cube_sum_given_cubes_given_a_e("../results/parameters.txt", "../results/cube_sum_vectors.txt");

	return 0;
}*/
