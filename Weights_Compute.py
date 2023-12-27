# This script aims to compute the coupling weights, namely, 

import numpy as np
from sympy.physics.wigner import wigner_3j

lmax = 5

def Odd(i,j,k):
    if (-1)**(i+j+k) == -1:
        return True
    else:
        return False

for i in range(lmax+1):
    for j in range(lmax+1):
        for k in range(abs(i-j), min(lmax+1, i+j+1)):
            if Odd(i,j,k):
                continue
            for m1 in range():
                for m2 in range():
                    for m3 in range():
                        





# Below are codes for debuging and testing
if __name__ == "__main__":
    1 == 1