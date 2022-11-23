# Practical cube-attack against nonce-misused ASCON
# Filename : CHK22_conditional_cubes_choices.py
# Date : May 2022
# Author : Jules Baudrin
# Content :  This file presents the sets mentioned in Figure 8 of our work. It helps understanding the differences
# between the choices made in [EPRINT:CHK22] and the choices we made regarding conditional cubes in the nonce-misuse
# setting.

from sage.all import *
from ascon_functions import *
import re


def track_multiplications_nonce_misuse(target_row, target_var, ring, ring2):
    # Track the multiplications occurring during the second S-box layer in the nonce-misuse setting and return all
    # the quadratic monomials present in the ANF for further studies

    # Initialize the state after the first S-box layer with only terms of degree 1
    Y = [[ring(0) for i in range(64)] for j in range(5)]
    for i in range(64):
        Y[0][i] = ring('v%d' % i) * ring('a%d + 1' % i)
        Y[1][i] = ring('v%d' % i)
        Y[3][i] = ring('v%d' % i) * ring('e%d' % i)
        Y[4][i] = ring('v%d' % i) * ring('a%d' % i)

    Y[0][0] = ring('v0')  # a_0 = 0
    Y[4][0] = 0

    Y = lin_layer(Y)  # first linear layer
    Y = sbox_layer(Y)  # second S-box layer

    # Regex in order to look for the targeted variables, either inside a term or at the end of a term
    last_var = re.compile('v%d$' % target_var)
    middle_var = re.compile('v%d\*' % target_var)

    # All quadratic v_0v_i as a dictionary mapping the index i to the list of coefficients of v_0v_i appearing in the ANF
    set_all_quad_terms = dict([(i, set([])) for i in range(64)])

    # for each coordinate, split the polynomial into a list of terms, recover only terms of degree 2 in public variables
    # which contains the primary variable
    for i in range(5):
        for j in range(64):
            coordinate_quad_terms = 0
            for term in str(Y[i][j]).split(' + '):
                if (middle_var.search(term) or last_var.search(term)) and (term.count('v') == 2):
                    coordinate_quad_terms += ring2(term)
            if coordinate_quad_terms:
                for (monom, coeff) in zip(coordinate_quad_terms.monomials(), coordinate_quad_terms.coefficients()):
                    set_all_quad_terms[int(str(monom).split('*')[1][1:])].add(coeff)

    return set_all_quad_terms


if __name__ == '__main__':
    r = BooleanPolynomialRing(192, ['v%d' % d for d in range(64)] +
                                 ['a%d' % d for d in range(64)] +
                                 ['e%d' % d for d in range(64)])

    r2 = PolynomialRing(PolynomialRing(GF(2), ['a%d' % d for d in range(64)] +
                                 ['e%d' % d for d in range(64)]), ['v%d' % d for d in range(64)])

    primary_variable = 0  # can be changed to i in {0, ..., 63} to shift the sets
    print('Primary variable considered:', primary_variable)

    print('Among the 64-1 = 63 other variables, WHEN a0 = 0:')
    indices_coeff_dict = track_multiplications_nonce_misuse(0, primary_variable, r, r2)
    nb_multiplied_variables = len([a for a in indices_coeff_dict if indices_coeff_dict[a]])
    print('\t- ', 63 - nb_multiplied_variables, 'are NEVER multiplied by v_0 after 2 rounds.')
    print('\t- ', nb_multiplied_variables, 'are multiplied by v_0 after 2 rounds.')

    print('\nAmong the %d variables multiplied by v_0:'%nb_multiplied_variables)
    var_without_cst_coeff = [a for a in indices_coeff_dict if indices_coeff_dict[a] and not(1 in indices_coeff_dict[a])]
    print('\t- ', nb_multiplied_variables - len(var_without_cst_coeff),
          'have at least one coefficient constant and equal to 1:')
    print('\tNo condition on the inner state will force the disappearance of all x_0x_i in this set.\n')
    print('\t- ', len(var_without_cst_coeff), 'DO NOT HAVE a coefficient constant and equal to 1:\n\t',
          var_without_cst_coeff)

    print('\nAmong those %d variables with non constant coefficients:'%len(var_without_cst_coeff))
    variables_with_non_cst_gcd = [a for a in var_without_cst_coeff if gcd(indices_coeff_dict[a]) != 1]
    variables_with_e0_gcd = [a for a in var_without_cst_coeff if r2('e0').divides(gcd(indices_coeff_dict[a]))]
    variables_with_non_e0_gcd = [a for a in variables_with_non_cst_gcd if not(a in variables_with_e0_gcd)]
    print('\t- %d have a GCD of the coefficients which is divisible by e_0 \n\t(adding e_0 = 0 to a_0 = 0 will force '
          'the disappear of all x_0x_i) :'%len(variables_with_e0_gcd), variables_with_e0_gcd)
    for a in variables_with_e0_gcd:
        print('\t\t', a, indices_coeff_dict[a], gcd(indices_coeff_dict[a]))

    print('\n\t- %d have GCD of the coefficients which is NOT divisible by e_0 but their GCDs are linear:\n\t'
          '(adding a linear condition (depending on i) to a_0 = 0 will force the disappear of all '
          'x_0x_i) :\n\t'%len(variables_with_non_e0_gcd), variables_with_non_e0_gcd)
    for a in variables_with_non_e0_gcd:
        print('\t\t', a, indices_coeff_dict[a], gcd(indices_coeff_dict[a]))

    var_with_complementary_coeff = []
    for a in var_without_cst_coeff:
        for coeff in indices_coeff_dict[a]:
            if (coeff + 1) in indices_coeff_dict[a]:
                var_with_complementary_coeff.append(a)
                break
    print('\n\t- %d have a constant GCD because at least one coefficient is the complementary of another (p and p+1) '
          '\n\t(No condition on the inner state will force the disappearance of all x_0x_i '
          'in this set.) :'%len(var_with_complementary_coeff), var_with_complementary_coeff)
    for a in var_with_complementary_coeff:
        print('\t\t', a, indices_coeff_dict[a])

    print('\n\t- For the remaining variables (26-7-12-2 = 5), more than one linear conditions on the inner state '
          'must be added\n\tto a_0 = 0 to force the disappearance of all x_0x_i')
    for a in var_without_cst_coeff:
        if not(a  in variables_with_e0_gcd or a in variables_with_non_e0_gcd or a in var_with_complementary_coeff):
            print('\t\t', a, indices_coeff_dict[a])