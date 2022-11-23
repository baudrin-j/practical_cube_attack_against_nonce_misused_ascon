# Practical cube-attack against nonce-misused ASCON
# Filename : conditional_cubes_choices.py
# Date : May 2022
# Author : Jules Baudrin
# Content : This file explains the choices made in our work regarding conditional cubes in the nonce-misuse setting.
# It shows how the sets S_{*} were built.

from sage.all import *
from ascon_functions import *
import re


def track_multiplications_nonce_misuse(target_row, primary_var):
    # Track the multiplications occurring during the second S-box layer in the nonce-misuse setting
    # primary_var is the variable tracked. primary_var coming from the target_row will have the flag T (for "target")
    # primary_var coming from other rows will have the flag O (for "other")
    
    ring = BooleanPolynomialRing(194, ['v%d' % d for d in range(64)] +
                                 ['a%d' % d for d in range(64)] +
                                 ['e%d' % d for d in range(64)] +
                                 ['T', 'O'])  # flags

    # Initialize the state after the first S-box layer with only terms of degree 1
    Y = [[ring(0) for i in range(64)] for j in range(5)]
    for i in range(64):
        Y[0][i] = ring('v%d' % i) * ring('a%d + 1' % i)
        Y[1][i] = ring('v%d' % i)
        Y[3][i] = ring('v%d' % i) * ring('e%d' % i)
        Y[4][i] = ring('v%d' % i) * ring('a%d' % i)

    # Add the flags to terms containing the primary_var
    for i in range(5):
        if i == target_row:
            Y[i][primary_var] *= ring('T')
        else:
            Y[i][primary_var] *= ring('O')

    Y = lin_layer(Y)  # first linear layer
    Y = sbox_layer(Y)  # second S-box layer

    s_tot = set([])  # s_tot will contain all variables multiplied by primary_var through the second S-box layer
    s_target = set([])  # s_target will contain all variables multiplied by some primary_var coming from target_row
    s_others = set([])  # s_others will contain all variables multiplied by some primary_var coming from another row

    # Regex in order to look for the targeted variables, either inside a term or at the end of a term
    last_var = re.compile('v%d$' % primary_var)
    middle_var = re.compile('v%d\*' % primary_var)

    # for each coordinate, split the polynomial into a list of terms, recover only terms of degree 2 in public variables
    # which contains the primary variable, and keep only the index of the 2nd variable.
    # Add each index to the proper set(s).
    for i in range(5):
        for j in range(64):
            for term in str(Y[i][j]).split(' + '):
                if (middle_var.search(term) or last_var.search(term)) and (term.count('v') == 2):
                    term_as_list = term.split('*')
                    if term_as_list[0] == ('v%d' % primary_var):
                        index = int(term_as_list[1][1:])
                    else:
                        index = int(term_as_list[0][1:])

                    s_tot.add(index)
                    if 'T' in term_as_list:
                        s_target.add(index)
                    else:
                        s_others.add(index)

    s_tot = set([i for i in range(64) if i != primary_var]).difference(s_tot)
    # s_tot now contains all variables not multiplied by the primary variable
    s_target = s_target.difference(s_others)
    # s_target now contains all variables multiplied by primary_variable from target row only

    return sorted(s_target), sorted(s_tot)


if __name__ == '__main__':
    primary_variable = 0  # can be changed to i in {0, ..., 63} to shift the sets.

    s = [[] for i in range(5)]
    s[0], _ = track_multiplications_nonce_misuse(0, primary_variable)
    s[2], _ = track_multiplications_nonce_misuse(3, primary_variable)
    s[1], s[3] = track_multiplications_nonce_misuse(4, primary_variable)
    # the second output of the three calls to the function is the same, we only keep the last one

    # From this knowledge we can easily compute S5:
    s[4] = sorted(set(range(64)).difference([0]).difference(s[0]).difference(s[1]).difference(s[2]).difference(s[3]))
    for i in range(5):
        print('S%d' % i, len(s[i]), s[i])
