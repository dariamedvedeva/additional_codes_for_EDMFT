from scipy.sparse import *
from scipy import *
import numpy as np
import os
import scipy.sparse.linalg as slin

def readfile(filename):
    ##############
    # read data
    ##############
    f = open(filename)
    matrix_file = f.readlines()

    if matrix_file.__len__() > 1:
        os.abort()

    # i - index
    i = 0
    j = 0

    array = []
    for element in matrix_file:
        value = element.split()
        length = value.__len__()
        for i in range(0, length, 2):
            array.append(np.complex(np.float(value[i]), np.float(value[i + 1])))
            j += 1

    f.close()
    return array

def readfile_col(filename):
    ##############
    # read indices
    ##############
    f = open(filename)
    matrix_file = f.readlines()

    if matrix_file.__len__() > 1:
        os.abort()

    # i - index
    i = 0
    j = 0

    array = []
    for element in matrix_file:
        value = element.split()
        length = value.__len__()
        for i in range(length):
            array.append(np.int(value[i]) - 1)
            j += 1
    f.close()
    return array

def readfile_row(filename):
    ##############
    # read index pointer
    ##############
    f = open(filename)
    matrix_file = f.readlines()

    # i - index
    i = 0
    j = 0

    array = []
    for element in matrix_file:
        value = element.split()
        length = value.__len__()
        # print "Number of elements in matrix in ", filename, length
        for i in range(length):
            array.append(np.int(value[i]) - 1)
            j += 1
    f.close()
    return array

def print_matrix_to_file(matrix, filename, print_matrix_in_file):
    if (print_matrix_in_file == 1):
        # print matrix to file
        with open(filename, 'w') as file:
            for i in range(numrows):
                for j in range(numrows):
                    row_in_matrix = matrix[i, j]
                    file.write(str(row_in_matrix))
                    file.write('\t')
                file.write('\n')
        file.close()

        filename_print_real_part = "real_" + filename
        with open(filename_print_real_part, 'w') as file:
            for i in range(numrows):
                for j in range(numrows):
                    row_in_matrix = np.round(real(matrix[i, j]), 2)
                    file.write(str(row_in_matrix))
                    file.write(' ')
                file.write('\n')
        file.close()

        filename_print_imag_part = "imag_" + filename
        with open(filename_print_imag_part, 'w') as file:
            for i in range(numrows):
                for j in range(numrows):
                    row_in_matrix = imag(matrix[i, j])
                    file.write(str(row_in_matrix))
                    file.write('\t')
                file.write('\n')
        file.close()
        print "Matrix was saved in files"


# data prepared
data = np.asarray(readfile("fort.444"))
col = np.asarray(readfile_col("fort.445"))
row_pntr = np.asarray(readfile_row("fort.446"))

# construction of the density matrix
A = csr_matrix((data, col, row_pntr))
matrix = A.todense()

# size of matrix
# need check of matrix shape
numrows = len(matrix)
print "Matrix dimention is ", numrows

# calculation of eigenpairs
vals, vecs = slin.eigsh(A, which='SA')
print "Eigen values:"
print vals

# ground state (minimum)
print "g.s. = ", matrix.min()

#conjugate matrix
conjugate_matrix = matrix.conj()
conjugate_matrix = conjugate_matrix.transpose()
# print conjugate_matrix

# check hermicity
for i in range(numrows):
    for j in range(numrows):
        if (matrix[i, j] != conjugate_matrix[i, j]):
            print 'trouble'
            os.abort()
    if(i == numrows - 1):
        print "Matrix is hermitian"

print_matrix_to_file(matrix, "matrix.dat", 1)