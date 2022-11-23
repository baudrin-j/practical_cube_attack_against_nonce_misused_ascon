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


uint value(uint64_t word, uint i) {
	return (word >> (63 - i)) & 1;
}


uint xor_values(vector<uint> values){
	bool x = 0;
	for(auto &y : values)
		x ^= y;
	return x;
}


/*
 * Given an input "parameters" file containing the already-recovered vectors a
 * and e, as well as a list of cubes of size 31 this functions :
 *  1- first selects pseudo-random values for vectors b and c;
 *  2- then, computes the cube-sum vectors corresponding to all cubes stored
 *  in the input file;
 *  3- and finally, outputs a file containing both the used values for b and c,
 *  as well as the cube-sum vectors.
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
	uint64_t b = random_monom();
	uint64_t c = random_monom();

	ofstream outputfile(outputfilename);
	outputfile << std::hex << b << endl << c << endl;
	for(uint i = 2; i < lines.size(); i++) {
		uint64_t cube_int = lines[i];
		vector<uint> cube;
		for(uint j = 0; j < 64; j++) {
			if((cube_int >> (63 - j)) & 1)
				cube.insert(cube.end(), j);
		}
		auto start = high_resolution_clock::now();
		uint64_t state[5] = {0,0,0,0,0};
		state[1] = a;
		state[2] = b;
		state[3] = c;
		state[4] = ~(c ^ e);

		cube_sum(state, rounds, cube, last_lin, cst);
		outputfile << std::hex << state[0] << endl;

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout << std::dec << "Time: " << duration.count() << endl;
	}
	outputfile.close();
}


int main(){

	omp_set_num_threads(8);
	cube_sum_given_cubes_given_a_e("../results/parameters.txt", "../results/cube_sum_vectors.txt");

	return 0;
}
