/*
 * Practical cube-attack against nonce-misused ASCON
 * Filename : permutation.cpp
 * Date : May 2022
 * Author : Jules Baudrin
 * Content : All the necessary functions to apply ASCON permutation to a 320-bit
 * state. A state is represented as an array of 5 uint64_t.
*/
#include "permutation.h"


// Shift function for 64-bit word
uint64_t rotr64(uint64_t x, unsigned int n)
{
	return (x >> n) | (x << (64 - n));
}


// Returns Sigma_i(x) (see ASCON specifications for more details)
uint64_t sigma(uint64_t x, unsigned int i)
{
	unsigned int a, b;
	switch(i)
   {
	case 0 :
		a = 19;
		b = 28;
		break;
	case 1 :
		a = 61;
		b = 39;
		break;
	case 2 :
		a = 1;
		b = 6;
		break;
	case 3 :
		a = 10;
		b = 17;
		break;
	case 4:
		a = 7;
		b = 41;
		break;

	default : return 0;
   }
	return x ^ rotr64(x, a) ^ rotr64(x, b);

}


// Bit-sliced ASCON sbox-layer, MODIFIES the state x
void sbox(uint64_t* x)
{
	uint64_t t0, t1, t2, t3, t4;
	x[0] ^= x[4];
	x[2] ^= x[1];
	x[4] ^= x[3];
	t0 = (~x[0]) & x[1];
	t1 = (~x[1]) & x[2];
	t2 = (~x[2]) & x[3];
	t3 = (~x[3]) & x[4];
	t4 = (~x[4]) & x[0];
	x[0] ^= t1;
	x[1] ^= t2;
	x[2] ^= t3;
	x[3] ^= t4;
	x[4] ^= t0;
	x[1] ^= x[0];
	x[3] ^= x[2];
	x[0] ^= x[4];
	x[2] = (~x[2]);
}


// ASCON constant addition, MODIFIES the state x
void add_cst(uint64_t* x, unsigned int i, unsigned int nb_rounds)
{
	x[2] ^=  ((i + 12 - nb_rounds) ^ ((15 - (i + 12 - nb_rounds)) << 4));
}


/*
 * ASCON permutation, MODIFIES the state x
 *   - Constant additions can be skipped or not depending on parameter cst.
 *   - Parameters i (current round index) and nb_rounds (total number of rounds)
 *     are given in order to use the appropriate constant sequence.
 *   - The last linear layer can be skipped or not depending on param lin_layer.
 */
void p(uint64_t* x, unsigned int i, unsigned int nb_rounds, bool lin_layer, bool cst)
{
	if(cst)
		add_cst(x, i, nb_rounds);
	sbox(x);
	if(lin_layer)
	{
		for(unsigned int j = 0; j < 5; j++)
			x[j] = sigma(x[j], j);
	}
}

/*
 * ASCON iterated permutation, MODIFIES the state x
 *   - Constant additions can be skipped or not depending on parameter cst.
 *   - nb_rounds is the total number of rounds: p^{nb_rounds} is applied to x.
 *   - The last linear layer can be skipped or not depending on param lin_layer.
 */
void multi_p(uint64_t* x, unsigned int nb_rounds, bool cst, bool last_linlayer)
{
	for(unsigned int i = 0 ; i < nb_rounds; i++)
		p(x, i, nb_rounds, (i != (nb_rounds - 1))  || last_linlayer, cst);
}
