from math import *
import numpy as np

#------Set the expected parameter here------#

H0 = 67.6 # Hubble constant
Omega_Mass = 0.31 # default Omega(matter)
Omega_Lambda = 0.0# default Omega(vacumm)


#------Functions to compute the distance from the redshift------#

def Redshift_to_Distance(z, type="DCMR", h = H0/100.0, Omega_m = Omega_Mass, n = 1000):
    Omega_Radiation = 4.165E-5/(h*h)
    Omega_total = Omega_m + Omega_Lambda + Omega_Radiation
    Omega_Curvature = 1 - Omega_total
    a = 1.0
    Tyr = 977.8
    c = 299792.458
    
    
    
    if(type == "DCMR"):
        
        az = 1.0/(1+1.0*z)
        age = 0.
        n=1000         # number of points in integrals
        for i in range(n):
            a = az*(i+0.5)/n
            adot = np.sqrt(Omega_Curvature+(Omega_m/a)+(Omega_Radiation/(a*a))+(Omega_Lambda*a*a))
            age = age + 1./adot

        zage = az*age/n
        zage_Gyr = (Tyr/H0)*zage
        DTT = 0.0
        DCMR = 0.0

        # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
        for i in range(n):
            a = az+(1-az)*(i+0.5)/n
            adot = np.sqrt(Omega_Curvature+(Omega_m/a)+(Omega_Radiation/(a*a))+(Omega_Lambda*a*a))
            DTT = DTT + 1./adot
            DCMR = DCMR + 1./(a*adot)

        DTT = (1.-az)*DTT/n
        DCMR = (1.-az)*DCMR/n
        age = DTT+zage
        age_Gyr = age*(Tyr/H0)
        DTT_Gyr = (Tyr/H0)*DTT
        DCMR_Gyr = (Tyr/H0)*DCMR
        DCMR_Mpc = (c/H0)*DCMR
        return(DCMR_Mpc)
        
    elif(type == "DA"):
        return(0)
      
    elif(type == "DL"):
        return(0)
        
    elif(type == ""):
        return(0)
        
    elif(type == ""):
        
        return(0)
    
    
# print(Redshift_to_Distance(0.5))