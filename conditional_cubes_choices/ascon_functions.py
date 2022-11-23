# Practical cube-attack against nonce-misused ASCON
# Filename : ascon_functions.py
# Date : May 2022
# Author : Jules Baudrin
# Content : This file contains the linear-layer and the sbox-layer functions used in the permutation of ASCON.


def lin_layer(X):
    A1 = X[0][-19:] + X[0][:-19]
    A2 = X[0][-28:] + X[0][:-28]

    B1 = X[1][-61:] + X[1][:-61]
    B2 = X[1][-39:] + X[1][:-39]

    C1 = X[2][-1:] + X[2][:-1]
    C2 = X[2][-6:] + X[2][:-6]

    D1 = X[3][-10:] + X[3][:-10]
    D2 = X[3][-17:] + X[3][:-17]

    E1 = X[4][-7:] + X[4][:-7]
    E2 = X[4][-41:] + X[4][:-41]

    for i in range(64):
        X[0][i] += A1[i] + A2[i]
        X[1][i] += B1[i] + B2[i]
        X[2][i] += C1[i] + C2[i]
        X[3][i] += D1[i] + D2[i]
        X[4][i] += E1[i] + E2[i]

    return X


def sbox_layer(X):
    for i in range(64):
        X[0][i] += X[4][i]
        X[4][i] += X[3][i]
        X[2][i] += X[1][i]

        x12 = (X[1][i] + 1) * X[2][i]
        x34 = (X[3][i] + 1) * X[4][i]
        x01 = (X[0][i] + 1) * X[1][i]
        x23 = (X[2][i] + 1) * X[3][i]
        X[3][i] += (X[4][i] + 1) * X[0][i]
        X[4][i] += x01
        X[0][i] += x12
        X[1][i] += x23
        X[2][i] += x34

        X[1][i] += X[0][i]
        X[3][i] += X[2][i]
        X[0][i] += X[4][i]
        X[2][i] += 1

    return X