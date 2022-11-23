# Practical cube-attack against nonce-misused ASCON
# Filename : system_solving.py
# Date : May 2022
# Author : Jules Baudrin
# Content : Sage script which uses a previously-created file "cube_sum_vectors.txt" in order to recover information
#           from the system solving. "cube_sum_vectors.txt" is output by "values_recovery".
#           Information recovered are verified and compared with the genuine values of b and c.

from sage.all import *
from sage.sat.boolean_polynomials import solve as solve_sat


def get_bin(byte, N):
    X = []
    for i in range(N):
        t = (byte >> (N-1-i)) & 1
        X.append(t)
    return X


if __name__ == '__main__':
    # Store the genuine values of b and c, as well as the cube-sum vectors of the targeted cube of size 31.
    f = open("results/cube_sum_vectors.txt")
    lines = []
    for line in f:
        print(line)
        lines.append(int("0x"+line, 16))
        print(lines[-1])
    b = lines[0]
    c = lines[1]
    b = get_bin(b, 64)
    c = get_bin(c, 64)

    ring = BooleanPolynomialRing(128, ['b%d'%d for d in range(64)] + ['c%d'%d for d in range(64)])
    dico = dict(zip(ring.gens(), b + c))

    polys = []
    error = 0
    count_cst_poly = 0

    # For each targeted cube:
    for i in range(2, len(lines)):
        # 1. Recover the cube-sum vector
        out = lines[i]
        out = get_bin(out, 64)

        # 2. Recover the computed coefficients
        f = open("results/polynomials_cube_%d.txt"%(i - 2))
        count = 0

        # For each output coefficient
        for poly in f:
            # 3. Build the polynomial, monomial by monomial
            poly = poly.split(' + ')
            p = ring(0)
            for m in poly:
                p += ring(m)

            # Check if p is constant, for statistics
            if p == ring(0) or p == ring(1):
                count_cst_poly += 1

            # 4. Add the corresponding cube value to build an equation
            p += ring(out[count])

            # Check if the recovered value matches the value of p(b, c). It should always be the case.
            if p.subs(dico):
                error += 1
                print(count, p)
            else:
                polys.append(p)

            count += 1
        f.close()

    print("\nNB error found:", error, "NB cst coefficients found:", count_cst_poly)

    # Solve the system (through cryptominisat)
    solution_sat = solve_sat(polys, n=infinity)

    print("Number of solutions for the current system:", len(solution_sat))  # Should be compared to 2^{128}

    # If a single solution is found, it should be the genuine pair of vectors (b, c)
    if len(solution_sat) == 1:
        b_recovered = [0 for i in range(64)]
        c_recovered = [0 for i in range(64)]
        for k,v in solution_sat[0].items():
            if str(k)[0] == 'b':
                b_recovered[int(str(k)[1:])] = v
            else:
                c_recovered[int(str(k)[1:])] = v
        print("Is the recovered b equal to the genuine b?: ", b == b_recovered)
        print("Is the recovered c equal to the genuine c?: ", c == c_recovered)
