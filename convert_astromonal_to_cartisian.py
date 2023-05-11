from astropy.io import fits
import numpy as np



# All in units of h-1 Mpc
H = 70
h = 0.7
c = 9.715611890751e-15 * h
c_normal = 299792.458

def Redshift_to_Distance(redshift):
    return h * (redshift * c_normal / H)


def Converts_Spherical_to_Cartesian(spherical):
    ra = spherical[:,2]
    dec = spherical[:,1]
    dis = spherical[:,0]
    
    x = dis * np.sin(dec) * np.cos(ra)
    y = dis * np.sin(dec) * np.sin(ra)
    z = dis * np.cos(dec)
    
    result = np.vstack((x,y,z)).T
    
    return result


def read_file(file_path):
    data_file = fits.open(file_path)
    
    ra = data_file[1].data.field('RA')
    dec = data_file[1].data.field('DEC')
    redshift = data_file[1].data.field('Z')
    
    distance = Redshift_to_Distance(redshift)
    
    spherical = np.vstack((distance, dec, ra)).T
    cartesian = Converts_Spherical_to_Cartesian(spherical)
    
    data_file.close()
    
    return cartesian
    
    
test_result = read_file("/Volumes/Liyang/Research/Data_Galaxy/galaxy_DR12v5_CMASS_South.fits")

print(test_result)

np.save('initial_data.npy', test_result)

