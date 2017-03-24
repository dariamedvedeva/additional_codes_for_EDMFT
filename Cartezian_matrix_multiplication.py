def test_example():
    # matrices should be square
    dimension1 = 3
    dimension2 = 2
    columns = dimension1
    rows = dimension1
    columns2 = dimension2
    rows2 = dimension2
    Matrix_1 = [[0 for x in range(rows)] for y in range(columns)]
    Matrix_2 = [[0 for x in range(rows2)] for y in range(columns2)]
    increment_element = 1
    for ii in range(rows):
        for jj in range(columns):
            Matrix_1[ii][jj] = increment_element
            increment_element += 1

    increment_element = 1
    for ii in range(rows2):
        for jj in range(columns2):
            Matrix_2[ii][jj] = increment_element + columns * rows
            increment_element += 1

    print_matrix("A = \n", Matrix_1)
    print_matrix("B = \n", Matrix_2)
    cartezian_matrix_myltipl(Matrix_1, Matrix_2)


def print_matrix(text, input):
    print text
    numrows = len(input)
    # numcols = len(input[0])
    for i in range(numrows):
        print input[i][:]
    print "\n"


def cartezian_matrix_myltipl(A, B):
    # intent -  A and B matrices
    # output - C
    # realization:
    # c_{i,j; k,l} = a_{i,k} * b_{j,l}
    numrows_A = len(A)
    numcols_A = len(A[0])
    i = numrows_A
    k = numcols_A
    numrows_B = len(B)
    numcols_B = len(B[0])
    j = numrows_B
    l = numcols_Bs
    C = [[0 for x in range(i * j)] for y in range(k * l)]
    # A
    for ii in range(i):
        for kk in range(k):
            # B
            n = 0
            for jj in range(j):
                m = 0
                for ll in range(l):
                    C[n + j * ii][m + l * kk] = A[ii][kk] * B[jj][ll]
                    m += 1
                n += 1
    print_matrix("C = \n", C)
    # In fact we obtain that second matrix is multiplied by every element
    # from the first matrix. And result matrix presents smthg like follow
    #
    #   a_11 * B      a_12 * B    a_13 * B    ...     a_1n * B
    #   a_21 * B      a_22 * B    a_23 * B    ...     a_2n * B
    #      ...          ...          ...      ...       ...
    #   a_21 * B      a_22 * B    a_23 * B    ...     a_2n * B

test_example()
