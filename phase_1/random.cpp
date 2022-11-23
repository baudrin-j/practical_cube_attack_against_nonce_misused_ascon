/*
 * Filename : random.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : Random uint64_t generation functions
*/
#include "random.h"

using namespace std;

/*
 * Generates and returns a random 64-bit word through the C++ pseudo-random
 * number generation library
 */
monom random_monom() {
	monom m = 0;
	std::random_device rd;  // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> distrib(0, 1);
	for(uint i = 0; i < 64; i++) {
		uint t = distrib(gen);
		if(t)
			m |= ((uint64_t) 1) << i;
	}
	return m;
}

/*
 * Generates and returns a random 64-bit word of fixed weight w through the C++
 * pseudo-random number generation library
 */
monom random_monom_weight(uint w) {
	monom m = 0;
	bool loop;
	uint t;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> distrib(0, 63);
	for(uint i = 0; i < w; i++) {
		loop = true;
		while(loop) {
			t = distrib(gen);
			if(((m >> t) & 1) == 0) {
				m |= (((uint64_t) 1) << t);
				loop = false;
			}
		}
	}
	return m;
}