# In this matrix, we want to import the edge correction algorithm to account for the survey geometry

import numpy as np
import sympy.physics.wigner as wg
import NPCF_estimator as Ne


l_max = 5

def D(l, L, l_prime):
    
    return ((2*l + 1)*(2 * L + 1)*(2 * l_prime + 1)) ** 0.5


def Coupling_Matrix(l1, l2, l3, l1_prime, l2_prime, l3_prime, bin_lists, Random_catalogs, weights_lists):
    
    sum = 0
    
    for i in range(-l_max, l_max + 1):
        for j in range(-l_max, l_max + 1):
            for k in range(-l_max, l_max + 1):
                sum += (Ne.Estimator(l1, l2, l3, Random_catalogs, weights_lists, bin_lists) / Ne.Estimator(0, 0, 0, Random_catalogs, weights_lists, bin_lists)) * D(l1, i, l1_prime) * D(l2, j, l2_prime) * D(l3, k, l3_prime) * wg.wigner_9j(l1, i, l1_prime, l2, j, l2_prime, l3, k, l3_prime, prec = 64)
                
    
    M = (4 * np.pi) ** -1.5 * sum
    
    return M


def Coupling_Matrix_generate(bin_lists, Random_catalogs, weights_lists):
    list = np.zeros(((l_max+1) ** 3, (l_max+1) ** 3))
    for i in range(l_max+1):
        for j in range(l_max+1):
            for k in range(l_max+1):
                row_index = i * l_max**2 + j * l_max + k
                
                for l in range(l_max+1):
                    for m in range(l_max+1):
                        for n in range(l_max+1):
                            column_index = l * l_max**2 + m * l_max + n

                            list[row_index, column_index] = Coupling_Matrix(i,j,k,l,m,n, bin_lists, Random_catalogs, weights_lists)
    return list

