# Practical cube-attack against nonce-misused ASCON
# Filename : system_solving.py
# Date : May 2022
# Author : Jules Baudrin
# Content : Sage script made to be launch from "coefficient_recovery" in order to recover information
#           from the system solving. Input and output files are stored in subfolder "results".
#           NB : PLEASE CHANGE THE VALUE OF VARIABLE "dir_results" BEFORE LAUNCHING THE PROGRAM.
#           Information recovered are verified and compared with the genuine values of a.

from sage.all import *
from sage.sat.boolean_polynomials import solve as solve_sat


def get_bin(byte, N):
    X = []
    for i in range(N):
        t = (byte >> (N-1-i)) & 1
        X.append(t)
    return X


if __name__ == '__main__':
    # NEED TO BE MODIFIED ACCORDINGLY
    # FOR EXAMPLE: dir_results = "/Users/usr/Documents/phase_2/results/"
    dir_results = "/Users/usr/Documents/phase_2/results/"

    # Recover the cube-sum vector value
    f = open(dir_results + "cube_sum_vectors.txt")
    cube_sum_vect = 0
    for line in f:
        print(line)
        cube_sum_vect = int("0x"+line, 16)
        cube_sum_vect = get_bin(cube_sum_vect, 64)
    f.close()

    # Recover the genuine vector a value
    f = open(dir_results + "parameters.txt")
    for line in f:
        a = int("0x"+line, 16)
        break
    a = get_bin(a, 64)
    f.close()

    ring = BooleanPolynomialRing(64, ['a%d'%d for d in range(64)])
    dico = dict(zip(ring.gens(), a))

    # Recover the polynomial expression of the coefficients
    polys = []
    error = 0
    count_cst_poly = 0
    f = open(dir_results + "polynomials.txt")
    count = 0
    for poly in f:
        poly = poly.split(' + ')
        p = ring(0)
        for m in poly:
            p += ring(m)

        # Verify if p is constant or not, for statistics
        if p == ring(0) or p == ring(1):
            count_cst_poly += 1

        # Add the cube-sum value to build the current equation
        p += ring(cube_sum_vect[count])
        polys.append(p)

        # Verify is the equation is correct or not, thanks to the genuine value of a
        # For verification only. No error message should be print.
        if p.subs(dico):
            error += 1
            print(count, p)
        else:
            polys.append(p)
        count += 1

    f.close()


    print("NB err:", error, "NB cst poly", count_cst_poly)
    # Solve the built system via CRYPTOMINISAT
    solution_sat = solve_sat(polys, n=infinity)

    # Checks if a unique solution exists. If it is the case, output the information to be used in the next iterations.
    if len(solution_sat) != 1:
        print("Underdefined system, information may be recovered with more tedious work.")
    else:
        f = open(dir_results + "recovered_a.txt", 'w')
        for k, v in solution_sat[0].items():
            index = int(str(k)[1:])
            print(a[index] == v)
            f.write(str(k) + " = " + str(v) + "\n")
        f.close()
