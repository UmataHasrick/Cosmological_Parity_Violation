import sys, os
import numpy as np
import multiprocessing
from sympy.physics.wigner import wigner_3j, wigner_9j
import subprocess
## First read-in the input parameters from the command line
#if len(sys.argv)!=5:
    #raise Exception("Need to specify N, LMAX, N_threads and all-parity flag")
#else:
N = int(4)
LMAX = int(3)
threads = int(10)#CPU相关的
all_parity = int(0) # if 1, also compute odd-parity weights
#0是只算even，1是even odd都算
print("\nComputing the edge-correction matrices for the %dPCF up to l_max = %d\n"%(N,LMAX))

# import timing module if present (no problem if not)


#给出l1,l2,l3的可能取值
if N==4:
    # First define array of ells
    ell_1,ell_2,ell_3 = [[] for _ in range(3)]#Three empty lists
    for l1 in range(0,LMAX+1,1):
        for l2 in range(0,LMAX+1,1):
            for l3 in range(abs(l1-l2),min(l1+l2,LMAX)+1,1):
                if not all_parity:
                    if (-1.)**(l1+l2+l3)==-1: continue#odd就不输出结果到数组，返回循环
                ell_1.append(l1)
                ell_2.append(l2)
                ell_3.append(l3)

    # Define coupling matrix, by iterating over all Lambda triples
    print("Computing 4PCF coupling matrix on %d CPUs"%threads)
def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def compute_matrix_coeff():
        # output submatrix
    Mcoefficient = np.zeros((len(ell_1),len(ell_1),len(ell_1)))
    for i in range(len(ell_1)):
        # i is first matrix index
        L_1,L_2,L_3=ell_1[i],ell_2[i],ell_3[i] # l1 l2 l3 which stored in the same index(i)of 3 different array
        pref_1 = np.sqrt((2.*L_1+1.)*(2.*L_2+1.)*(2.*L_3+1.))

        for j in range(len(ell_1)):
            # j is second matrix index
            Lpp_1,Lpp_2,Lpp_3=ell_1[j],ell_2[j],ell_3[j]# l1'' l2'' l3''(L1,L2,L3 in the paper) which stored in the same index(i)of 3 different array
            pref_2 = pref_1*np.sqrt((2.*Lpp_1+1.)*(2.*Lpp_2+1.)*(2.*Lpp_3+1.))
            # add phase
            pref_2 *= (-1.)**(Lpp_1+Lpp_2+Lpp_3)

            for k in range(len(ell_1)):
                # k indexes inner Lambda' term
                Lp_1,Lp_2,Lp_3 = ell_1[k],ell_2[k],ell_3[k]# l1' l2' l3' which stored in the same index(i)of 3 different array

                # Compute prefactor
                pref = pref_2*np.sqrt((2.*Lp_1+1.)*(2.*Lp_2+1.)*(2.*Lp_3+1.))/(4.*np.pi)**(3./2.)

                # Compute 3j couplings
                three_j_piece = np.float64(wigner_3j(L_1,Lp_1,Lpp_1,0,0,0))
                if three_j_piece==0: continue
                three_j_piece *= np.float64(wigner_3j(L_2,Lp_2,Lpp_2,0,0,0))
                if three_j_piece==0: continue
                three_j_piece *= np.float(wigner_3j(L_3,Lp_3,Lpp_3,0,0,0))
                if three_j_piece==0: continue

                # Compute the 9j component
                nine_j_piece = np.float64(wigner_9j(L_1,Lp_1,Lpp_1,L_2,Lp_2,Lpp_2,L_3,Lp_3,Lpp_3,prec=8))
                if nine_j_piece==0: continue

                Mcoefficient[i,j,k] =  pref  * three_j_piece * nine_j_piece#这里没求和，只是算了系数出来，最后还是要乘上f_L1,L2,L3
    return Mcoefficient
Mcoefficient=compute_matrix_coeff()
print(Mcoefficient)
a=Mcoefficient
with open('my_file.txt', 'w') as f:
    # Convert the list to a string and write it to the file
    f.write(str(a))
print(a)




# 这就算了个grant给我