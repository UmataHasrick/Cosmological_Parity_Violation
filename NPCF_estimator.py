import numpy as np
import wigner as wg
import scipy.special as spe
import math as m


def longitude(x,y):
    
    """Return the longitude of a given point (x,y,z) within range [0,2*PI). The two args must be of the same length. If not, the program will exit with return value 0. 

    Args:
        x (NDArray): One-dimensional array of the x component(s) of given point(s).
        y (NDArray): One-dimensional array of the y component(s) of given point(s)

    Returns:
        ndarray: One-dimensional array of the longitude of the points within range [0,2*PI)
    """

    output = np.empty(0)

    length_x = x.shape[0]
    length_y = y.shape[0]

    assert length_x == length_y

    for i in range(length_x):
        value = m.atan2(y[i],x[i])
        
        if value < 0:
            value += 2 * np.pi

        output = np.append(output, value)

    return output


def Converts_cartesian_to_spherical(r_cartesian):
    
    """The function can convert the list of the position(s) of the point(s) in cartesian form into spherical form with in range [0,Infinity), [0,PI) and [0,2*PI)
    
    Args:
        r_cartesian (NDArray): (-1,3) array of the position(s) of the point(s) in cartesian form. 

    Returns:
        NDArray: (-1,3) array of the position(s) of the point(s) in cartesian form with in range [0,Infinity), [0,PI) and [0,2*PI) (A.k.a in the form of (r, theta, phi))
    """

    x = r_cartesian[:,0]
    y = r_cartesian[:,1]
    z = r_cartesian[:,2]

    r = (x ** 2 + y ** 2 + z **2) ** 0.5 

    theta = np.arccos(z/r)
    phi = longitude(x,y)

    return np.transpose(np.nan_to_num(np.array([r, theta, phi])))


def bin_func(x, bin_min, bin_max):
    
    """The function can serve as a bin function to distinguish whether the input is within the bin.
    
    Args:
        x (Float): The input parameter to be distinguished. 
        b_min (Float): The lower limit of the bin. 
        b_max (Float): The upper limit of the bin

    Returns:
        Int: Return 1 if the x is with in the range and return 0 otherwise (Of the C Bool style)
    """
    
    if (x < bin_max) and (x > bin_min):
        return 1
    else:
        return 0
    

def sphere_volumn(r):
    
    """The function can calculate the volumn of a sphere of given radius. 
    
    Args:
        r (Float or NDArray): The radius(radii) of the sphere(s) to be calculated. 

    Returns:
        Float or NDArray: The volumn of the sphere(s) to be calculated
    """
    
    return (4 * np.pi * r ** 3) / 3


def V_b(bin_min, bin_max):
    
    """The function of the volumn of the spherical shell. It is used as the denominator to calculated the bin average. 

    Args:
        bin_min (Float or NDArray): The lower limit(s) of the bin(s). 
        bin_max (Float or NDArray): The upper limit(s) of the bin(s)
    
    Returns:
        Float: The volumn(s) of the spherical shell(s)
    """
    
    return sphere_volumn(bin_max) - sphere_volumn(bin_min)


def wigner_3j(l1, l2, l3, m1, m2, m3):
    
    """_summary_

    Returns:
        _type_: _description_
    """
    
    if m1 + m2 + m3 != 0:
        return 0
    else:
        intermedia = wg.wigner_3jm(l1, l2, l3, m1)
        m2_min = intermedia[0]
        m2_max = intermedia[1]
        value_list = intermedia[2]

        index = int(m2 - m2_min)

        return value_list[index]
    

def C_Epsilon_Lambda(l1, l2, l3, m1, m2, m3):
    return wigner_3j(l1, l2, l3, m1, m2, m3)


def A_func(l, m, primary_vertices, secondary_vertices, bin, weights_lists, average_caculate = True, configuration_space_average = False):

    Num_of_vertices = secondary_vertices.shape[0]
    sum = 0
    v_b = 1
    b_min = bin[0]
    b_max = bin[1]

    relative_position = secondary_vertices - primary_vertices
    relative_position_spherical = Converts_cartesian_to_spherical(relative_position)

    # print(relative_position_spherical)

    if (average_caculate == True) and (configuration_space_average == True):

        for j in range(Num_of_vertices): 
            if relative_position_spherical[j][0] < b_max and relative_position_spherical[j][0] > b_min:
                sum += weights_lists[j] * spe.sph_harm(m, l, relative_position_spherical[j][2], relative_position_spherical[j][1])
                # print(sum)

        v_b = v_b * V_b(b_min, b_max)

        return sum / v_b
    
    elif (average_caculate == True) and (configuration_space_average == False):
        
        for j in range(Num_of_vertices): 
            if relative_position_spherical[j][0] < b_max and relative_position_spherical[j][0] > b_min:
                sum += relative_position_spherical[j][0] * weights_lists[j] * spe.sph_harm(m, l, relative_position_spherical[j][2], relative_position_spherical[j][1])
                # print(sum)

        v_b = v_b * V_b(b_min, b_max)

        return sum / v_b
        
    elif (average_caculate == False) and (configuration_space_average == False):
        
        for j in range(Num_of_vertices): 
            if relative_position_spherical[j][0] < b_max and relative_position_spherical[j][0] > b_min:
                sum += relative_position_spherical[j][0] * weights_lists[j] * spe.sph_harm(m, l, relative_position_spherical[j][2], relative_position_spherical[j][1])
                # print(sum)
                
        v_b = V_b(b_min, b_max)

        return sum
        
    else:
    
        for j in range(Num_of_vertices): 
            if relative_position_spherical[j][0] < b_max and relative_position_spherical[j][0] > b_min:
                sum += weights_lists[j] * spe.sph_harm(m, l, relative_position_spherical[j][2], relative_position_spherical[j][1])
                # print(sum)

        return sum


def cache_IO(cache_file, information, num, IO_status, flag = 'a'):
    
    if IO_status == 0:
        IO_file = open(cache_file,flag)
        print(str(num)+"\t"+information, file = IO_file)
    
    
def Estimator(l1, l2, l3, Data_catalog, Data_catalog_weights_lists, bin_lists, Random_catalog, Random_catalog_weights_lists, DmR_status, space_volume,  cache_flag = 0, average_calculate = True, configuration_space_average = False):
    
    
    sum = 0


    if DmR_status == 0:
        # DDDD
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary,configuration_space_average=configuration_space_average, average_caculate=average_calculate) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary, configuration_space_average=configuration_space_average, average_caculate=average_calculate) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary, configuration_space_average=configuration_space_average, average_caculate=average_calculate)
                        sum += sum_terms
                    
    elif DmR_status == 1:
        # DDDR
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, Random_catalog, bin_lists[2], Random_catalog_weights_lists)
                        sum += sum_terms

    elif DmR_status == 2:
        # DDRD
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, Random_catalog, bin_lists[1], Random_catalog_weights_lists) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
    
    elif DmR_status == 3:
        # DRDD
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Random_catalog, bin_lists[0], Random_catalog_weights_lists) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
                        
    elif DmR_status == 4:
        # RDDD
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Data_catalog, bin_lists[0], Data_catalog_weights_lists) * A_func(l2, j, primary, Data_catalog, bin_lists[1], Data_catalog_weights_lists) * A_func(l3, k, primary, Data_catalog, bin_lists[2], Data_catalog_weights_lists)
                        sum += sum_terms
                        
    elif DmR_status == 5:
        # DDRR
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, Random_catalog, bin_lists[1], Random_catalog_weights_lists) * A_func(l3, k, primary, Random_catalog, bin_lists[2], Random_catalog_weights_lists)
                        sum += sum_terms
                        
                        
    elif DmR_status == 6:
        # DRRD
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Random_catalog, bin_lists[0], Random_catalog_weights_lists) * A_func(l2, j, primary, Random_catalog, bin_lists[1], Random_catalog_weights_lists) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
                        
                        
    elif DmR_status == 7:
        # DRDR
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Random_catalog, bin_lists[0], Random_catalog_weights_lists) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, Random_catalog, bin_lists[2], Random_catalog_weights_lists)
                        sum += sum_terms
                        
                        
    elif DmR_status == 8:
        # RRDD
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, Data_catalog, bin_lists[1], Data_catalog_weights_lists) * A_func(l3, k, primary, Data_catalog, bin_lists[2], Data_catalog_weights_lists)
                        sum += sum_terms
                        
    elif DmR_status == 9:
        # RDRD
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Data_catalog, bin_lists[0], Data_catalog_weights_lists) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, Data_catalog, bin_lists[2], Data_catalog_weights_lists)
                        sum += sum_terms
    
    elif DmR_status == 10:
        # RDDR
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, Data_catalog, bin_lists[1], Data_catalog_weights_lists) * A_func(l3, k, primary, Data_catalog, bin_lists[2], Data_catalog_weights_lists)
                        sum += sum_terms
    
    elif DmR_status == 11:
        # DRRR
        Num_of_Data_catalog = Data_catalog.shape[0]
        Num_of_weights = Data_catalog_weights_lists.shape[0]
        
        assert Num_of_Data_catalog == Num_of_weights

        for l in range(Num_of_Data_catalog):
            primary = Data_catalog[l]
            secondary = np.delete(Data_catalog, l, axis = 0)
            weights_primary = Data_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Data_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Random_catalog, bin_lists[0], Random_catalog_weights_lists) * A_func(l2, j, primary, Random_catalog, bin_lists[1], Random_catalog_weights_lists) * A_func(l3, k, primary, Random_catalog, bin_lists[2], Random_catalog_weights_lists)
                        sum += sum_terms
                        
    elif DmR_status == 12:
        # RDRR
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, Data_catalog, bin_lists[0], Data_catalog_weights_lists) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
                        
    elif DmR_status == 13:
        # RRDR
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, Data_catalog, bin_lists[1], Data_catalog_weights_lists) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
                        
    elif DmR_status == 14:
        # RRRD
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, Data_catalog, bin_lists[2], Data_catalog_weights_lists)
                        sum += sum_terms
    
    elif DmR_status == 15:
        # RRRR
        Num_of_Random_catalog = Random_catalog.shape[0]
        Num_of_weights = Random_catalog_weights_lists.shape[0]
        
        assert Num_of_Random_catalog == Num_of_weights

        for l in range(Num_of_Random_catalog):
            primary = Random_catalog[l]
            secondary = np.delete(Random_catalog, l, axis = 0)
            weights_primary = Random_catalog_weights_lists[l]
            weights_lists_secondary = np.delete(Random_catalog_weights_lists, l, axis = 0)
            
            for i in range(-l1, l1+1):
                for j in range(-l2, l2+1):
                    for k in range(-l3, l3+1):
                        sum_terms = weights_primary * C_Epsilon_Lambda(l1, l2, l3, i, j, k) * A_func(l1, i, primary, secondary, bin_lists[0], weights_lists_secondary) * A_func(l2, j, primary, secondary, bin_lists[1], weights_lists_secondary) * A_func(l3, k, primary, secondary, bin_lists[2], weights_lists_secondary)
                        sum += sum_terms
               
    return(sum/space_volume)

