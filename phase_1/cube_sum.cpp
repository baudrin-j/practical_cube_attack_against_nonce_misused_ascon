/*
 * Filename : cube_sum.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : parallelized cube-sum computation
*/
#include "cube_sum.h"
using namespace std;

/*
 * Computes the cube sum of a given cube.
 * - cube_index is the list of the cube variables.
 * - rounds is the number of rounds.
 * - partial_init is the given initial state. Only the four last rows
 * (inner state) matter.
 * - partial_init IS MODIFIED to return the cube_sum
 */
void cube_sum(uint64_t* partial_init, const uint &rounds, \
		const vector<uint> &cube_index, const bool &last_linlayer, \
		const bool &cst)
{
	uint64_t sum0 = 0;
	uint64_t sum1 = 0;
	uint64_t sum2 = 0;
	uint64_t sum3 = 0;
	uint64_t sum4 = 0;

	// subset is incidence vector of any subset of [0, cube_index.size() - 1]
#pragma omp parallel for default(none) shared(partial_init, rounds, last_linlayer, cst, std::cout, cube_index) reduction(^: sum0)  reduction(^: sum1)  reduction(^: sum2)  reduction(^: sum3)  reduction(^: sum4)
	for(uint64_t subset = 0; subset < (((uint64_t) 1) << cube_index.size()); subset++) {
		uint64_t state[5] = {0, 0, 0, 0, 0};

		// subset_mask is the corresponding subset of the cube
		uint64_t subset_mask = 0;
		for(uint i = 0; i < 64; i++) {
			if((subset >> i) & 1)
				subset_mask |= (((uint64_t) 1) << (63 - cube_index[i]));
		}

		// 1.a Initialize with the current vector of public variables.
		state[0] = subset_mask;
		// 1.b Use the given inner-state
		for(uint i = 1; i < 5; i++) { // a, b, c, d
			state[i] = partial_init[i];
		}

		// 2. Compute x rounds of ASCON
		multi_p(state, rounds, cst, last_linlayer);

		// 3. Sum up.
		sum0 ^= state[0];
		sum1 ^= state[1];
		sum2 ^= state[2];
		sum3 ^= state[3];
		sum4 ^= state[4];
	}

	// Copy the cube sum in partial_init
	partial_init[0] = sum0;
	partial_init[1] = sum1;
	partial_init[2] = sum2;
	partial_init[3] = sum3;
	partial_init[4] = sum4;
}
