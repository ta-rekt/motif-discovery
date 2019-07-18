# L-8 MCS 507 Fri 14 Sep 2012 : permanent.py

"""
Given an n-by-n matrix A, the permanent of the matrix is
the sum of the products A[i,s[i]], for i = 1,2,..,n,
where the sum runs over all permutations s of (1,2,..,n).
The expansion formula for the permanent above is very similar
to the row expansion formula for the determinant, except for
the sign changes, which are absent in the permanent expansion.
The Python function permanent below takes on input
a matrix and that returns the permanent of the matrix.
An auxiliary recursive function uses the expansion formula
above to compute the permanent.
"""

import numpy as np

def per(mtx, column, selected, prod, output=False):
    """
    Row expansion for the permanent of matrix mtx.
    The counter column is the current column,
    selected is a list of indices of selected rows,
    and prod accumulates the current product.
    """
    if column == mtx.shape[1]:
        if output:
            print (selected, prod)
        return prod
    else:
        result = 0
        for row in range(mtx.shape[0]):
            if not row in selected:
                result = result \
                + per(mtx, column+1, selected+[row], prod*mtx[row,column])
        return result

def permanent(mat):
    """
    Returns the permanent of the matrix mat.
    """
    return per(mat, 0, [], 1)
    

def main():
    """
    Test on the permanent.
    """
    dim = int(input('give the dimension : '))
    rmt = np.random.random_integers(0, 1, size=(dim, dim))
    print ('a random 0/1-matrix :\n', rmt)
    print ('permanent :', permanent(rmt))

if __name__ == "__main__":
    main()
