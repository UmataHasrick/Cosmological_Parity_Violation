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
LMAX = int(2)
threads = int(10)#CPU相关的
all_parity = int(1) # if 1, also compute odd-parity weights
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
    tmp_out = np.zeros((len(ell_1),len(ell_1),len(ell_1)))
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

                tmp_out[i,j,k] =  pref  * three_j_piece * nine_j_piece#其实这里没求和，只是算了系数出来，最后还是要乘上f_L1,L2,L3
    return tmp_out

Cmatrix=compute_matrix_coeff()
inputs = str(sys.argv[1])
##########################################



R_file = inputs+'.r_%dpcf.txt'%N
if N>2:#for different radial bins, we have different R and N
    countsR_original_data = np.loadtxt(R_file,skiprows=5+N) # skipping rows with radial bins,不包含radial bins 那几行

    # Extract ells and radial bins
if N==4:
    ell_1,ell_2,ell_3 = np.asarray(countsR_original_data[:,:3],dtype=int).T#前三列，就是l1,l2,l3三个量子数
    max_ell = np.max(ell_1)
    countsR = countsR_original_data[:,3:]#第四列以后的，2D的矩阵
    bin1,bin2,bin3 = np.loadtxt(R_file,skiprows=6,max_rows=3)#跳过前六行读取下面三行
    # Now load in D-R pieces and average
countsN_all = []
total_DmR = 0

for i in range(100):
    DmR_file = inputs+'.n%s_%dpcf.txt'%(str(i).zfill(2),N)##
    if not os.path.exists(DmR_file): continue
        # Extract counts
        # skipping rows with radial bins and ell
    if N==4:
        countsN_all.append(np.loadtxt(DmR_file,skiprows=9)[:,3:])#跟前面一样跳过前9行和前三列，3D的矩阵，后面跳列的操作好像是错的,这里得改，根据Nfile的内容
       
countsN_all = np.asarray(countsN_all)
#N_files = len(countsN_all)
countsN = np.mean(countsN_all,axis=0)#get mean of the quantum numbers of each bin [mean_bin1,mean_bin2,mean_bin3,....]
#countsN=countsN_all



if N==4:
        # Define coupling coefficients, rescaling by R_{Lambda=0}
    assert ell_1[0]==ell_2[0]==ell_3[0]#调试验证
    f_Lambda = countsR/countsR[0] # (first row should be unity!)#f_L1,L2,L3 的算法

    # Decide if we're using all-parity or even-parity multiplets
    if(np.sum((-1)**(ell_1+ell_2+ell_3)==1)==len(ell_1)):
            # no odd-parity contributions!
        all_parity=0
    else:
            # contains odd-parity terms!
        all_parity=1

        # Load coupling matrix
    LMAX = max(ell_1)
    if all_parity:
        input_weights = get_script_path()+'/../coupling_matrices/edge_correction_matrix_%dpcf_LMAX%d_all.npy'%(N,LMAX)
    else:
        input_weights = get_script_path()+'/../coupling_matrices/edge_correction_matrix_%dpcf_LMAX%d.npy'%(N,LMAX)
    if os.path.exists(input_weights):
        print("Loading edge correction weights from file.")
    else:
            # Compute weights from scratch
        subprocess.run(["python",get_script_path()+"/edge_correction_weights.py","%d"%N,"%d"%LMAX,'%d'%threads,'%d'%all_parity])

    coupling_tmp = Cmatrix

        # Define coupling matrix, by iterating over all Lambda triples
    coupling_matrix = np.zeros((len(ell_1),len(ell_1),len(bin1)))
    for i in range(len(ell_1)):
        for j in range(len(ell_1)):
            for k in range(len(ell_1)):
                coupling_matrix[i,j] += coupling_tmp[i,j,k]*f_Lambda[k]#这里把L求和了 生成真正的M,生成一整行的[sum,sum,sum.......]

        ## Now invert matrix equation to get zeta
        # Note that our matrix definition is symmetric
    zeta = np.zeros_like(countsN)
    for i in range(len(bin1)):
        zeta[:,i] = np.matmul(np.linalg.inv(coupling_matrix[:,:,i]),countsN[:,i]/countsR[0,i])#计算4PCF-estimator，利用矩阵乘法算M^-1乘上N_l'/R_0
print(zeta)

        # Now save the output to file, copying the first few lines from the N files





# 这就算了个grant给我
