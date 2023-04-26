# In this matrix, we want to import the edge correction algorithm to account for the survey geometry

import numpy as np
import sympy.physics.wigner as wg
import NPCF_estimator as Ne


l_max = 5

def D(l, L, l_prime):
    
    return ((2*l + 1)*(2 * L + 1)*(2 * l_prime + 1)) ** 0.5


def Coupling_Matrix(l1, l2, l3, l1_prime, l2_prime, l3_prime, r1, r2, r3, bin, R):
    
    sum = 0
    
    for i in range(-l_max, l_max):
        for j in range(-l_max, l_max):
            for k in range(-l_max, l_max):
                sum += Ne.Estimator(l1, l2, l3, R, )
                
    
    M = (4 * np.pi) ** -1.5
    
    return M