import numpy as np
import random

# The programme aims to provide a solution to generate the random positions for at least 1,500 primary points. 
# The distribution of the points are rather uniform. 

def Num_1D_pick(point_list, Num, lower_limit, upper_limit, minimal_sep):

    Num -= 1

    if Num >= 0:
        standby = random.uniform(lower_limit, upper_limit)
        point_list.append(standby)

        distance = upper_limit - lower_limit
        distance_left = standby - lower_limit
        distance_right = upper_limit - standby

        # print("\n")
        # print("The lower limit is %f, the upper limit is %f. " % (lower_limit, upper_limit))
        # print("We choose %f. " % (standby))
        # print("\n")

        if distance_left <= minimal_sep and distance_right <= minimal_sep:
            return Num
        
        elif distance_left <= minimal_sep and distance_right > minimal_sep:
            return Num_1D_pick(point_list, Num, standby + minimal_sep, upper_limit, minimal_sep)
    
        elif distance_right <= minimal_sep and distance_left > minimal_sep:
            return Num_1D_pick(point_list, Num, lower_limit, standby - minimal_sep, minimal_sep)
        else:
            Num_left = int(Num * ((distance_left - minimal_sep) / (distance - 2 * minimal_sep)))
            Num_right = Num - Num_left
            m1 = Num_1D_pick(point_list, Num_left, standby + minimal_sep, upper_limit, minimal_sep)
            m2 = Num_1D_pick(point_list, Num_right, lower_limit, standby - minimal_sep, minimal_sep)
            return m1 + m2
    
    else:
        return 0
    
def Distance(r1, r2):
    dr = r1 - r2
    return np.linalg.norm(dr)

def Minimal_sep_test(test, list, minimal_sep):
    # This function is used to test whether the test is within or without the minimal_sep range of every object in the list. 
    l = (list.shape)[0]

    result = True

    for i in range(l):
        if Distance(test, list[i]) < minimal_sep:
            result = False
            break
    
    return result


def Num_3D_pick(point_list, Num, space, minimal_sep):
    space_x_lower = space[0][0]
    space_x_upper = space[0][1]
    space_y_lower = space[1][0]
    space_y_upper = space[1][1]
    space_z_lower = space[2][0]
    space_z_upper = space[2][1]

    while(Num > 0):
        x = random.uniform(space_x_lower, space_x_upper)
        y = random.uniform(space_y_lower, space_y_upper)
        z = random.uniform(space_z_lower, space_z_upper)

        standby = np.array([x,y,z])
        # print("standby in (%f, %f, %f)" % (x,y,z))
        if Minimal_sep_test(standby, point_list, minimal_sep):
            point_list = np.append(point_list, standby.reshape((1,3)), axis=0)
            Num -= 1
    
    return point_list






# point_list = []
# list = Num_1D_pick(point_list, 1500, 0, 1, 60 / (1000 * 3 ** 0.5))

# print(list)
# print(point_list)

# point_list = np.empty((0,3))
# space = np.array([[0,1000],[0,1000],[0,1000]])
# Num = 1500
# minimal_sep = 60

# point_list = Num_3D_pick(point_list, Num, space, minimal_sep)
# print(point_list.shape)