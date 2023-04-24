import matplotlib.pyplot as plt
import numpy as np
import random


def Minimal_sep_test(test, list, minimal_sep):
    # This function is used to test whether the test is without the minimal_sep range of every object in the list. 
    l = (list.shape)[0]

    result = True

    for i in range(l):
        if abs(test-list[i]) < minimal_sep:
            result = False
            break
    
    return result


def Non_repeating_random_numbers(Num_of_numbers, Minimal_sep_x):
    number = Num_of_numbers
    result_set = np.array([random.random()])

    while number > 0:
        pre_random_num = random.random()
        if Minimal_sep_test(pre_random_num, result_set, Minimal_sep_x):
            result_set = np.append(result_set, pre_random_num)
            number -= 1
            print(pre_random_num)
    return result_set


def Random_Pick_vertices(Space, Num_of_vertices, Minimal_sep):
    # Setting the parameters
    length_x = Space[0]
    length_y = Space[1]
    length_z = Space[2]
    Minimal_sep_normalized_x = Minimal_sep / (length_x * 3**0.5)
    Minimal_sep_normalized_y = Minimal_sep / (length_y * 3**0.5)
    Minimal_sep_normalized_z = Minimal_sep / (length_z * 3**0.5)

    print(Minimal_sep_normalized_x)

    # Generating the random numbers of the vertices
    lx = length_x * Non_repeating_random_numbers(Num_of_vertices, Minimal_sep_normalized_x)
    ly = length_y * Non_repeating_random_numbers(Num_of_vertices, Minimal_sep_normalized_y)
    lz = length_z * Non_repeating_random_numbers(Num_of_vertices, Minimal_sep_normalized_z)

    output = np.transpose(np.array([lx, ly, lz]))

    return  output 

Random_Pick_vertices([1000,1000,1000], 100, 40)