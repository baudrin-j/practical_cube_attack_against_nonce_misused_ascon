# Practical cube-attack against nonce-misused ASCON
# Filename : LDW17_conditional_cubes_choices.py
# Date : May 2022
# Content : This file explains the choices made in [ToSC:LDW17] regarding conditional cubes for 5-round and 6-round
# initializations. It uses the ANF point of view and the same reasoning that we introduced in our paper.
# It also enables to build more cubes with the same properties.

from sage.all import *
from ascon_functions import *
import re

# Tuples containing, the cube indices and the variables input on each nonce row given in [ToSC:LDW17].
# We also add a targeted row, depending on the ANF of the primary column after the first S-box layer

# INITIALIZATION : primary variable on row 4, all other variables on row 4 as well
# ANF OF THE PRIMARY COLUMN AFTER FIRST S-BOX LAYER:
#       k0 * v0 + IV * k0 + k0 * k64 + IV + k0 + k64
#       v0 + k0 * k64 + IV + k0 + k64
#       v0 + k0 + k64 + 1
#       (IV + 1) * v0 + IV + k0 + k64
#       (k0 + 1) * v0 + IV * k0 + k0

# targeted row = 0, so k0 is the corresponding coefficient
cube1_deg16 = ([0, 5, 8, 12, 14, 15, 18, 19, 21, 27, 28, 30, 34, 37, 49, 50], [], list(range(64)), 0)

# targeted row = 4, so (k0 + 1) is the corresponding coefficient
cube2_deg16 = ([0, 5, 7, 8, 14, 15, 24, 27, 30, 34, 37, 41, 43, 49, 50, 52], [], list(range(64)), 4)

#######################################################################################################################
# INITIALIZATION : primary variable on row 3, all other variables on row 3 as well
# NB: for cube5_deg32, not enough variables are available with this setting, that is why some other variables
# are input on row 4  as well by LDW. This choice limits the diffusion of those variables during the first round.
# ANF OF THE PRIMARY COLUMN AFTER FIRST S-BOX LAYER:
#       v0 + IV * k0 + k0 * k64 + IV + k0 + k64
#       (k0 + k64 + 1) * v0 + k0 * k64 + IV + k0 + k64
#       k0 + k64 + 1
#       (IV + 1) * v0 + IV + k0 + k64
#       v0 + IV * k0 + k0

# targeted row = 1, so (k0 + k64 + 1) is the corresponding coefficient
cube4_deg16 = ([0, 5, 8, 14, 15, 16, 17, 20, 27, 29, 30, 33, 34, 35, 37, 38], list(range(64)), [], 1)
cube4_deg32 = ([0, 1, 3, 4, 5, 6, 8, 14, 15, 16, 17, 20, 26, 27, 29, 30,
                33, 34, 35, 37, 38, 39, 40, 46, 49, 50, 55, 58, 59, 60, 62, 63], list(range(64)), [], 1)

# targeted row = 1, so (k0 + k64 + 1) is the corresponding coefficient
cube5_deg16 = ([0, 5, 8, 14, 15, 27, 29, 30, 34, 36, 37, 38, 39, 45, 49, 50], list(range(64)), [], 1)
cube5_deg32 = [[0, 1, 4, 5, 6, 8, 9, 11, 14, 15, 16, 18, 24, 26, 27, 29, 30, 34,
                36, 37, 38, 39, 45, 47, 48, 49, 50, 56, 58, 59, 60, 63], list(range(64)), [9, 11, 18, 24, 47], 1]

#######################################################################################################################
# INITIALIZATION : row 3 = public var, row 4 = null, except for the column of the primary variable
# NB: for cube1_deg32 and cube2_deg32 inputting more variables on row 4 may not be necessary ?
# ANF OF THE PRIMARY COLUMN AFTER FIRST S-BOX LAYER:
#       (k0 + 1) * v0 + IV * k0 + k0 * k64 + IV + k0 + k64
#       (k0 + k64) * v0 + k0 * k64 + IV + k0 + k64
#       k0 + k64 + 1
#       IV + k0 + k64
#       k0 * v0 + IV * k0 + k0

# targeted row = 1, so (k0 + k64) is the corresponding coefficient
cube3_deg16 = ([0, 1, 4, 5, 6, 8, 14, 15, 16, 17, 20, 26, 27, 29, 30, 33], list(range(64)), [0], 1)
cube3_deg32 = ([0, 1, 4, 5, 6, 8, 14, 15, 16, 17, 20, 26, 27, 29, 30, 33,
                34, 35, 37, 38, 39, 40, 46, 48, 49, 50, 55, 56, 58, 59, 62, 63], list(range(64)), [0], 1)

# targeted row = 4, so k0 is the corresponding coefficient
cube1_deg32 = ([0, 1, 4, 5, 6, 7, 8, 10, 13, 14, 15, 16, 17, 24, 26, 27,
                30, 34, 35, 37, 40, 41, 43, 46, 48, 49, 50, 52, 56, 59, 60, 63], list(range(64)), [0, 60, 63], 4)

# targeted row = 0, so k0 + 1 is the corresponding coefficient
cube2_deg32 = ([0, 1, 4, 5, 6, 8, 9, 12, 14, 15, 16, 17, 18, 19, 21, 26,
                27, 28, 30, 34, 35, 37, 40, 46, 48, 49, 50, 53, 56, 59, 60, 63], list(range(64)), [0, 60, 63], 0)
#######################################################################################################################


def get_bin(byte, N):
    X = []
    for i in range(N):
        t = (byte >> (N-1-i)) & 1
        X.append(t)
    return X


def track_multiplications_init(input_vars_row_3, input_vars_row_4, iv, target_row, prim_var, print_bool):
    # Track the multiplications occurring during the second S-box layer for ASCON initialization.
    #       - prim_var is the primary variable.
    #       - target_row is the row of origin of the primary variable.
    #       - input_vars_row_3 is the list of indices of the input variables on row 3.
    #       - input_vars_row_4 is the list of indices of the input variables on row 4.
    #       - print_bool is a printing option
    # prim_var coming from the target row will have the flag T (for "target")
    # target var coming from other rows will have the flag O (for "other")

    ring = BooleanPolynomialRing(194, ['v%d' % d for d in range(64)] + ['k%d' % d for d in range(128)] + ['T', 'O'])

    # Initialize the state as in ASCON initialization
    Y = [[ring(0) for i in range(64)] for j in range(5)]
    for i in range(64):
        Y[0][i] = iv[i]
        Y[1][i] = ring('k%d' % i)
        Y[2][i] = ring('k%d' % (i + 64))
        if i in input_vars_row_3:
            Y[3][i] = ring('v%d' % i)
        if i in input_vars_row_4:
            Y[4][i] = ring('v%d' % i)

    Y = sbox_layer(Y)  # first S-box layer

    # Add the flags to terms containing the primary variable
    for i in range(5):
        if i == target_row:
            dico = dict(zip([ring('v%d' % prim_var)], [ring('v%d*T' % prim_var)]))
            Y[i][prim_var] = Y[i][prim_var].subs(dico)
        else:
            dico = dict(zip([ring('v%d' % prim_var)], [ring('v%d*O' % prim_var)]))
            Y[i][prim_var] = Y[i][prim_var].subs(dico)

    Y = lin_layer(Y)  # first linear layer
    Y = sbox_layer(Y)  # second S-box layer

    s_tot = set([])  # s_tot will contain all variables multiplied by the primary variable after the second round
    s_target = set([])  # s_T will contain all variables multiplied by some primary variable coming from target_row
    s_others = set([])  # s_O will contain all variables multiplied by some v0 coming from another row than target_row

    # Regex in order to look for the targeted variables, either inside a term or at the end of a term
    last_var = re.compile('v%d$' % prim_var)
    middle_var = re.compile('v%d\*' % prim_var)

    # for each coordinate, split the polynomial into a list of terms, recover only terms of degree 2 in public variables
    # which contains the primary variable, and keep only the index of the 2nd variable.
    # Add each index to the proper set(s).
    for i in range(5):
        for j in range(64):
            t = 0
            for term in str(Y[i][j]).split(' + '):
                if (middle_var.search(term) or last_var.search(term)) and (term.count('v') == 2):
                    term_as_list = term.split('*')  # a term viewed as a list of variables
                    if term_as_list[0] == ('v%d' % prim_var):
                        index = int(term_as_list[1][1:])
                    else:
                        index = int(term_as_list[0][1:])

                    s_tot.add(index)
                    if 'T' in term_as_list:
                        s_target.add(index)
                    else:
                        s_others.add(index)

    s_tot = set([i for i in range(64) if i != prim_var]).difference(s_tot)
    # s_tot now contains all variables *NOT* multiplied by any primary variable

    s_target = s_target.difference(s_others)
    # s_target now contains all variables multiplied some primary variable coming from the target row only

    if(print_bool):
        print('Variables v_i only multiplied by some primary variable v_%d coming from row %d:\n'%(prim_var, target_row),
              len(s_target), sorted(s_target))
        print('\nVariables v_i that are never multiplied by v_%d:\n'%(prim_var),  len(s_tot), sorted(s_tot))
    return sorted(s_target), sorted(s_tot)


if __name__ == '__main__':
    choice_of_cube = cube1_deg16  # Can take the value cube_{1, 2, 3, 4, 5}_deg{16, 32}
    iv_vect = get_bin(0x80400c0600000000, 64)  # Can also take the value [0 for i in range(64)]

    # Using a null IV enables to *over-estimate* the diffusion of each variable through the first S-box layer:
    #       - with the genuine IV, if iv_i = 1, then variable v_i is absent from the fourth row,
    #       - with the all-zero IV iv_i = 0 for all i and thus v_i is always present on the fourth row.
    # With a null IV we obtain smaller sets s_target and s_tot: less conditional cubes can be found but the ones found
    # can be shifted cyclically. This is not true for all conditional cubes while targeting the initialization.
    # So a null IV can help to understand better the choice of cubes.

    print('Cube : ', choice_of_cube[0], '\n')
    if choice_of_cube == cube5_deg16 or choice_of_cube == cube5_deg32:  # Necessary shift of variables because cube5
        primary_var = 1                                                 # is never used for primary variable 1
        cube = [(a + primary_var) % 64 for a in choice_of_cube[0]]      # (see Algorithm 1 in [ToSC:LDW17])
        if choice_of_cube == cube5_deg32:
            choice_of_cube[2] = [(a + primary_var) % 64 for a in choice_of_cube[2]]
    else:
        primary_var = 0
        cube = choice_of_cube[0]

    s_target_row, s_not_mult = track_multiplications_init(choice_of_cube[1], choice_of_cube[2], iv_vect,
                                                          choice_of_cube[3], primary_var, True)

    intersection_target_row = set(s_target_row).intersection(cube)
    intersection_not_multiplied = set(s_not_mult).intersection(cube)
    print('\nChoices of variables among the two sets:')
    print(len(intersection_target_row), sorted(intersection_target_row))
    print(len(intersection_not_multiplied), sorted(intersection_not_multiplied))