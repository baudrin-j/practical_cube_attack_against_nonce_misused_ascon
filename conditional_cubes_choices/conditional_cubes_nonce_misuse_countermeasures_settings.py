# Practical cube-attack against nonce-misused ASCON
# Filename : conditional_cubes_counter_measures.py
# Date : May 2022
# Author : Jules Baudrin
# Content :  This file presents the different settings proposed as possible counter-measures in our paper.
# For each of the four settings, it helps to verify if the same logic used in conditional_cubes_choices.py can be used
# to build conditional cubes.

from sage.all import *
from ascon_functions import *
import re


def track_multiplications_nonce_misuse(target_row, primary_var, row_ordering):
    # Track the multiplications occurring during the second S-box layer in the nonce-misuse setting
    # primary_var is the tracked primary variable.
    #   - primary_var coming from the target_row will have the flag T (for "target")
    #   - primary_var coming from other rows will have the flag O (for "other")
    # row_ordering is the list of labels of input variables: for instance if row_ordering = ['a', 'v', 'b', 'c', 'd']
    # then we consider that the public variables are input on the second row

    ring = BooleanPolynomialRing(322, ['v%d' % d for d in range(64)] +
                                 ['a%d' % d for d in range(64)] +
                                 ['b%d' % d for d in range(64)] +
                                 ['c%d' % d for d in range(64)] +
                                 ['d%d' % d for d in range(64)] +
                                 ['T', 'O'])

    # Initialize the state depending on the ordering of variables given as parameter
    Y = [[ring(0) for i in range(64)] for j in range(5)]
    for i in range(64):
        for j in range(5):
            Y[j][i] = ring(row_ordering[j]+'%d'%i)

    Y = sbox_layer(Y)  # first S-box layer

    # Add the flags to terms containing the primary_var
    for i in range(5):
        if i == target_row:
            dico = dict(zip([ring('v%d' % primary_var)], [ring('v%d*T' % primary_var)]))
            Y[i][primary_var] = Y[i][primary_var].subs(dico)
        else:
            dico = dict(zip([ring('v%d' % primary_var)], [ring('v%d*O' % primary_var)]))
            Y[i][primary_var] = Y[i][primary_var].subs(dico)

    Y = lin_layer(Y)  # first linear layer
    Y = sbox_layer(Y)  # second S-box layer

    s_tot = set([])  # s_tot will contain all variables multiplied by v0 through the second S-box layer
    s_target = set([])  # s_T will contain all variables multiplied by some v0 coming from target_row
    s_others = set([])  # s_O will contain all variables multiplied by some v0 coming from another row than target_row

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
    # s_tot now contains all variables not multiplied by v0
    s_target = s_target.difference(s_others) # s_T now contains all variables multiplied by v0 from target row only

    return sorted(s_target), sorted(s_tot)


def anf_col(x):
    # ASCON S-box
    y = []
    y += [x[4] * x[1] + x[3] + x[2] * x[1] + x[2] + x[1] * x[0] + x[1] + x[0]]
    y += [x[4] + x[3] * x[2] + x[3] * x[1] + x[3] + x[2] * x[1] + x[2] + x[1] + x[0]]
    y += [x[4] * x[3] + x[4] + x[2] + x[1] + 1]
    y += [x[4] * x[0] + x[4] + x[3] * x[0] + x[3] + x[2] + x[1] + x[0]]
    y += [x[4] * x[1] + x[4] + x[3] + x[1] * x[0] + x[1]]
    return y


def print_anf_half_round(pos):
    # Print the ANF of a column after the first S-box layer,
    # depending on the index "pos" of the input row for public vars.
    R = BooleanPolynomialRing(5, ['v', 'a', 'b', 'c', 'd'])
    A = PolynomialRing(GF(2), ['a', 'b', 'c', 'd'])
    B = PolynomialRing(A, 'v')
    init_col = [R('a'), R('b'), R('c'), R('d')]
    init_col = init_col[:pos] + [R('v')] + init_col[pos:]
    col = anf_col(init_col)
    for l in col:
        print(B(str(l)))


def print_sets(row_ordering):
    # Print the interesting sets of variables to build conditional cubes,
    # depending on the chosen input row for public vars.

    primary_variable = 0  # can be changed to i in {0, ..., 63} to shift the sets
    s = [[] for i in range(6)]
    s[0], _ = track_multiplications_nonce_misuse(0, primary_variable, row_ordering)
    s[1], _ = track_multiplications_nonce_misuse(1, primary_variable, row_ordering)
    s[2], _ = track_multiplications_nonce_misuse(2, primary_variable, row_ordering)
    s[3], _ = track_multiplications_nonce_misuse(3, primary_variable, row_ordering)
    s[4], s[5] = track_multiplications_nonce_misuse(4, primary_variable, row_ordering)
    # the second output of the five calls to the function is the same, we only keep the last one

    for i in range(6):
        if i in range(5):
            print("Variables multiplied by some v_%d coming only from row %d: " % (primary_variable, i), end="")
        else:
            print("Variables never multiplied by any v%d: " % primary_variable, end="")
        print(len(s[i]), s[i])

        def print_anf_half_round(pos):
            # Print the ANF of a column after the first S-box layer,
            # depending on the index "pos" of the input row for public vars.
            R = BooleanPolynomialRing(5, ['v', 'a', 'b', 'c', 'd'])
            A = PolynomialRing(GF(2), ['a', 'b', 'c', 'd'])
            B = PolynomialRing(A, 'v')
            init_col = [R('a'), R('b'), R('c'), R('d')]
            init_col = init_col[:pos] + [R('v')] + init_col[pos:]
            col = anf_col(init_col)
            for l in col:
                print(B(str(l)))


if __name__ == '__main__':
    for i in range(5):
        order = ['a', 'b', 'c', 'd']
        order = order[:i] + ['v'] + order[i:]
        print('----------------------------------------')
        print('Current order:', order)
        print('----------------------------------------')
        print_anf_half_round(i)
        print_sets(order)
        print('\n')
