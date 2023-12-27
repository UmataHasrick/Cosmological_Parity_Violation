import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
from Redshift_and_Distance import Redshift_to_Distance

def plot_spherical_coordinates(theta, phi, r):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # 创建Basemap对象
    m = Basemap(projection='moll', lat_0=0, lon_0=0, resolution='l')

    # 将球坐标转换为经纬度坐标
    lon = phi
    lat = theta

    # 绘制球坐标图
    x, y = m(lon, lat)
    m.scatter(x, y, c=r, s = 5, cmap='viridis', edgecolors='black', linewidths=0.5)
    
    # 绘制参考线
    m.drawmeridians(range(-180, 180, 30), labels=[False, False, False, True])
    m.drawparallels(range(-90, 90, 30), labels=[True, False, False, False])
    
    # 绘制图例
    cbar = plt.colorbar()
    cbar.set_label('Distance')


    plt.show()

# 示例数据
# theta = np.linspace(0, np.pi, 50)  # theta 角度范围从 0 到 pi
# phi = np.linspace(0, 2*np.pi, 50)  # phi 角度范围从 0 到 2*pi
# theta, phi = np.meshgrid(theta, phi)  # 创建 theta 和 phi 的网格
# r = np.sin(theta)  # 示例中的 r 值为 sin(theta)

df = pd.read_csv("/Volumes/Liyang/Research/Data_Galaxy/Data_Decompressed/galaxy_DR12v5_CMASS_North.csv")

r = Redshift_to_Distance(df['Z'])
theta = df['DEC']
phi = df['RA']



plot_spherical_coordinates(theta, phi, r)