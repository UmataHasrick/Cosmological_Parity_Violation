import numpy as np
import matplotlib.pyplot as plt
import Figure_plotter
import Non_Repeating_random_numbers_generator as NRrng
import Tetrahedron_generator
import NPCF_estimator as Ne
import scipy.special as spe
import pandas as pd
import Redshift_and_Distance as RD
import matplotlib as mpl
from tqdm.auto import tqdm
import multiprocessing
import os
import logging
import sys
import datetime



# Set cosmological constant here
h = 0.676


input_data_dir = 
output_result_dir = 
cache_dir = 
file_list = 
log_dir = 

Num_of_threads = 

cache_IO_flag = False

lmax = 5
bin_schemes = [[],
               [],
               []]
bin_dict = {}

# For logging of the programme. Not important
# 创建日志文件名
log_file = log_dir+'A_computation_{}.log'.format(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

# 配置日志记录器
logger = logging.getLogger('A_computation')
logger.setLevel(logging.DEBUG)

# 创建文件处理器
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.DEBUG)

# 创建控制台处理器
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

# 创建日志格式器
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

# 设置文件处理器的格式
file_handler.setFormatter(formatter)

# 设置控制台处理器的格式
console_handler.setFormatter(formatter)

# 将处理器添加到日志记录器中
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# 记录日志
# logger.debug('This is a debug message')
# logger.info('This is an info message')
# logger.warning('This is a warning message')
# logger.error('This is an error message')
# logger.critical('This is a critical message')

for file in file_list:
    # Convert data file to proper format
    logger.info("File " + file + " now is getting xyzw. ")
    
    if cache_IO_flag == False or os.path.exists(cache_dir+file+"fancy.csv"):
        logger.info("File " + file + " converting data now")
        data = pd.read_csv(input_data_dir+file+".csv")
        
        Distance = RD.Redshift_to_Distance(data["z"].values) * h
        RA = data["RA"].values * np.pi / 180
        DEC = data["DEC"].values * np.pi / 180
        
        x = Distance * np.cos(RA) * np.cos(DEC)
        y = Distance * np.sin(RA) * np.cos(DEC)
        z = Distance * np.sin(DEC)
        
        W_FKP = data['WEIGHT_FKP'].values
        W_SYS = data['WEIGHT_SYSTOT'].values
        W_NOZ = data['WEIGHT_NOZ'].values
        W_CP = data['WEIGHT_CP'].values
        
        weights = W_FKP * W_SYS * (W_NOZ + W_CP - 1)
        
        pd.DataFrame(np.vstack((x,y,z,weights)).T,columns = ['X','Y','Z','W']).to_csv(cache_dir+file+"fancy.csv")
        
        del DEC,RA,W_FKP,W_CP,W_NOZ,W_SYS
    
    else:
        logger.info("Using cached data. ")
        data = pd.read_csv(cache_dir+file+"fancy.csv")
        x = data['X'].values
        y = data['Y'].values
        z = data['Z'].values
        weights = data['W'].values
    
    vertices = np.vstack((x,y,z)).T    
    del data,x,y,z
    
    # Now computing all the as'
    
    logger.info("Computing all the as' of {0} given the xyzw data. ".format(file))
    logger.info("lmax set to {0}. ".format(lmax))
    
    results = []
    
    for i in range(len(bin_schemes)):
        bin = bin_schemes[i]
        
        for l in range(0,lmax+1):
            for m in range(0,l+1):
                logger.info("Computing as' for l = {0}, m = {1}, bin scheme {2} - {3}. ".format(l,m,bin[0],bin[1]))
                logger.info("Establishing directory.")
                try:
                    os.mkdir(output_result_dir + "/" + file)
                except:
                    logger.info("Directory already established.")
                logger.info("Directory established. ")
                pool = multiprocessing.Pool(processes=Num_of_threads)
                
                with open(output_result_dir + "/" + file+"/" + file + "lm="+str(l)+str(m)+"_"+"bin="+str(bin[0])+"-"+str(bin[1])+".txt", mode='w') as w:
                    logger.info("File established. ")
                    for k in range(len(vertices)):
                        primary_vertices = vertices[k]
                        secondary_vertices = np.delete(vertices, k, axis = 0)
                        secondary_weights = np.delete(weights, k, axis = 0)
                        result = pool.apply_async(Ne.A_func, args=(l,m,primary_vertices, secondary_vertices, bin, secondary_weights), kwds={'average_caculate':True, 'configuration_space_average':False})
                        results.append(result)
                        
                    pool.close()
                    logger.info("Computing as' on {0} threads".format(Num_of_threads))
                    pool.join()
                    
                    for result in results:
                        w.write(str(np.imag(result.get())))
                        w.write("\n")
                
                logger.info("Computation as' for l = {0}, m = {1}, bin scheme {2} - {3} completed. ".format(l,m,bin[0],bin[1]))




file_handler.close()    