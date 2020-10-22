# try out ryser's formula for computing the permanent

import numpy as np
import itertools as it

# assumes input is square matrix
def permanent_fast(mat):

    indices = [i for i in range(mat.shape[0])]
    perm = 0
    sign = (-1)**len(indices)

    # loop numbers of columns (k)
    for k in range(0, mat.shape[1]):

        Ek = 0

        for c in it.combinations(indices, k):
            select = [i for i in indices if i not in c]

            Ak = np.take(mat, select, 1)

            PAk = np.multiply.reduce(np.add.reduce(Ak, 1))
            Ek += ((-1)**len(select))*PAk

        perm += Ek

    perm = perm*sign

    return perm
