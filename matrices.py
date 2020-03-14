import numpy as np

def MatxVec(MatA, Vecb, Dimd):
    C = []
    num = 0
    for i in range(Dimd):
        for j in range(Dimd):
            num += MatA[i, j] * Vecb[j]
        C.append(num)
    C = np.asarray(C)
    C = np.transpose(C)
    return C

def restaVec(Veca, Vecb, Dimd):
    c = []
    for i in range(Dimd):
        c.append(Veca[i] - Vecb[i])
    c = np.asarray(c)
    return c
