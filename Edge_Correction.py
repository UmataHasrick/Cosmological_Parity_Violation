# In this matrix, we want to import the edge correction algorithm to account for the survey geometry

import numpy as np
import sympy.physics.wigner as wg
import NPCF_estimator as Ne
from sympy.physics.wigner import wigner_3j, wigner_9j


def list_gen(odd_mode = True, lmax = 5):
    
    ell_1 = np.empty(0, dtype = np.int8)
    ell_2 = np.empty(0, dtype = np.int8)
    ell_3 = np.empty(0, dtype = np.int8)

    for i in range(0,lmax+1):
        for j in range(0,lmax+1):
            for k in range(abs(i-j),min(i+j,lmax)+1,1):
                
                if odd_mode == True:
                    if (-1.)**(i+j+k)==1: continue
                ell_1 = np.append(ell_1,i)
                ell_2 = np.append(ell_2,j)
                ell_3 = np.append(ell_3,k)
                
    return np.hstack((ell_1.reshape(-1,1), ell_2.reshape(-1,1), ell_3.reshape(-1,1)))

                    

def D(l, L, l_prime):
    
    return ((2.*l + 1.)*(2. * L + 1.)*(2. * l_prime + 1.)) ** 0.5



def Epsilon_Graunt_compute(Lambda_list):
    
    Num_of_list = Lambda_list.shape[0]
    EGc = np.zeros((Num_of_list, Num_of_list, Num_of_list))
    
    for i in range(Num_of_list):
        l1, l2, l3 = Lambda_list[i,1], Lambda_list[i,2], Lambda_list[i,3]
        coef_1 = D(l1, l2, l3)
        for j in range(Num_of_list):
            lpp1, lpp2, lpp3 = Lambda_list[j,1], Lambda_list[j,2], Lambda_list[j,3]
            coef_2 = -1.**(lpp1+lpp2+lpp3) * D(lpp1, lpp2, lpp3)
            for k in range(Num_of_list):
                lp1, lp2, lp3 = Lambda_list[k,1], Lambda_list[k,2], Lambda_list[k,3]
                coef_3 = D(lp1, lp2, lp3)
                
                EGc[i,j,k] = coef_1 * coef_2 * coef_3 * np.float64(wigner_3j(l1,lp1,lpp1,0,0,0)) * np.float64(wigner_3j(l2,lp2,lpp2,0,0,0)) * np.float64(wigner_3j(l3,lp3,lpp3,0,0,0)) * np.float64(wigner_9j(l1,lp1,lpp1,l2,lp2,lpp2,l3,lp3,lpp3,prec=8))

    return EGc

def f_func(bin_list, Random_catalog_weights, Random_catalog, space_volume, lambda_list):
    Num_of_list = lambda_list.shape[0]
    ff = np.zeros(Num_of_list)
    R0 = Ne.Estimator(0,0,0,Random_catalog, Random_catalog_weights, bin_list, np.zeros(1), np.zeros(1), 0, space_volume)
    for i in range(Num_of_list):
        Ri = Ne.Estimator(lambda_list[i,1], lambda_list[i,2], lambda_list[i,3], Random_catalog, Random_catalog_weights, bin_list, np.zeros(1), np.zeros(1), 0, space_volume)
        ff[i] = Ri/R0
    return ff
    

def Coupling_Matrix(f_func, EGc):
    
    assert f_func.shape[0] == EGc.shape[0]
    
    Num_of_list = EGc.shape[0]
    CM = np.zeros((Num_of_list, Num_of_list))
    
    for i in range(Num_of_list):
        for j in range(Num_of_list):
            for k in range(Num_of_list):
                CM[i,j] += f_func[k] * EGc[i,j,k]
                
    return CM
    
    
    
    
    
    
    
#     sum = 0
    
#     for i in range(-l_max, l_max + 1):
#         for j in range(-l_max, l_max + 1):
#             for k in range(-l_max, l_max + 1):
#                 sum += (Ne.Estimator(l1, l2, l3, Random_catalogs, weights_lists, bin_lists) / Ne.Estimator(0, 0, 0, Random_catalogs, weights_lists, bin_lists)) * D(l1, i, l1_prime) * D(l2, j, l2_prime) * D(l3, k, l3_prime) * wg.wigner_9j(l1, i, l1_prime, l2, j, l2_prime, l3, k, l3_prime, prec = 64)
                
    
#     M = (4 * np.pi) ** -1.5 * sum
    
#     return M


# def Coupling_Matrix_generate(bin_lists, Random_catalogs, weights_lists):
#     list = np.zeros(((l_max+1) ** 3, (l_max+1) ** 3))
#     for i in range(l_max+1):
#         for j in range(l_max+1):
#             for k in range(l_max+1):
#                 row_index = i * l_max**2 + j * l_max + k
                
#                 for l in range(l_max+1):
#                     for m in range(l_max+1):
#                         for n in range(l_max+1):
#                             column_index = l * l_max**2 + m * l_max + n

#                             list[row_index, column_index] = Coupling_Matrix(i,j,k,l,m,n, bin_lists, Random_catalogs, weights_lists)
#     return list

