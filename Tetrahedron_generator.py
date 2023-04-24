import numpy as np
import random

def Create_Single_Tetrahedron(position, parity, r, deviation):
    # In this function, it can build a single tetrahedron given the position of the primary point, the distance between the 
    # secondary points to the primary points, the parity of the tetrahedron, and the random deviation. 

    # position: a 3-d array, denoting the xyz value. 
    # parity: a two-value parameter only can take the value 1 and -1. 1 denotes the counterclockwise-constructed tetrahedron while the -1 denotes the clockwised ones. 
    # r: a 3-d array, denoting the 3 distance between the 3 secondary points to the primary points. The distances are listed from the shortest to the longest. 
    # deviation: a 6-d array. The first 3 value denotes the rotate angle of the vertex from the xyz axis and the last 3 values denote the deviation of the distance. 


    r_1 = r[0]
    r_2 = r[1]
    r_3 = r[2]
    Delta_angle_1 = deviation[0]
    Delta_angle_2 = deviation[1]
    Delta_angle_3 = deviation[2]
    Delta_r1 = deviation[3]
    Delta_r2 = deviation[4]
    Delta_r3 = deviation[5]

    primary = position

    transform_matrix_x = np.array([[        1,                          0,                      0],
                                   [        0,      np.cos(Delta_angle_1),  np.sin(Delta_angle_1)],
                                   [        0,     -np.sin(Delta_angle_1), np.cos(Delta_angle_1)]])
    transform_matrix_y = np.array([[np.cos(Delta_angle_2), 0, np.sin(Delta_angle_2)],
                                   [0,1,0],
                                   [-np.sin(Delta_angle_2), 0, np.cos(Delta_angle_2)]])
    transform_matrix_z = np.array([[np.cos(Delta_angle_3),  np.sin(Delta_angle_3),0],
                                   [-np.sin(Delta_angle_3), np.cos(Delta_angle_3),0],
                                   [0,0,1]])
    transform_matrix = transform_matrix_x @ transform_matrix_y @ transform_matrix_z

    r1_before_trans = parity * np.array([r_1 + Delta_r1,0,0])
    r2_before_trans = np.array([0,r_2 + Delta_r2,0])
    r3_before_trans = np.array([0,0,r_3 + Delta_r3])

    secondary_1 = primary + transform_matrix @ r1_before_trans
    secondary_2 = primary + transform_matrix @ r2_before_trans
    secondary_3 = primary + transform_matrix @ r3_before_trans

    return np.vstack((primary, secondary_1, secondary_2, secondary_3))




def Generate_random_deviations(r1_min, r1_max, r2_min, r2_max, r3_min, r3_max):
    # To generate a random deviation within the range. 

    Delta_angle_1 = random.uniform(-np.pi, np.pi)
    Delta_angle_2 = random.uniform(-np.pi, np.pi)
    Delta_angle_3 = random.uniform(-np.pi, np.pi)
    Delta_r1 = random.uniform(r1_min, r1_max)
    Delta_r2 = random.uniform(r2_min, r2_max)
    Delta_r3 = random.uniform(r3_min, r3_max)

    return np.array([Delta_angle_1, Delta_angle_2, Delta_angle_3, Delta_r1, Delta_r2, Delta_r3])




def Create_Multiple_Tetrahedron(vertices, parity, r, deviation_range):
    Num_of_vertices = (vertices.shape)[0]
    delta_r1_min = deviation_range[0]
    delta_r1_max = deviation_range[1]
    delta_r2_min = deviation_range[2]
    delta_r2_max = deviation_range[3]
    delta_r3_min = deviation_range[4]
    delta_r3_max = deviation_range[5]

    multiple_tetrahedron = np.zeros(shape = (0,4,3))

    for i in range(Num_of_vertices):
        deviation = Generate_random_deviations(delta_r1_min, delta_r1_max, delta_r2_min, delta_r2_max, delta_r3_min, delta_r3_max)
        position = vertices[i]
        single_tetrahedron = Create_Single_Tetrahedron(position, parity, r, deviation)
        multiple_tetrahedron = np.append(multiple_tetrahedron, single_tetrahedron.reshape((1,4,3)), axis = 0)

    return multiple_tetrahedron