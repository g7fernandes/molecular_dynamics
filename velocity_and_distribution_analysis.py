import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd 
import os
from glob import glob
import configparser

dirname = os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

config = configparser.ConfigParser()
config.read(res_dir + 'settings.txt')
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])
n_files = int(config['out_files']['out_files'].split()[0])

subspaces = float(input("Enter in how many subspaces the region of calculus will be divided.\nThe \
subspaces should contain particles enough for the statistic.\n"))
print("The velocity and distribution of the particles will be calculated in one direction.\n")
#direction = input("Enter the direction x or y\n")
print('x-direction\n')
direction = 'x'

step = input('Choose step:\n')
while int(step) >= 0: 
    pos = pd.read_csv(res_dir+"position.csv."+step, header=None, names = ["x","y"])
    vel = pd.read_csv(res_dir+"velocity.csv."+step, header=None, names = ["v_x", "v_y"])
    n = [x for x in range(len(pos))]
    n = pd.DataFrame(n, columns=["n"]) # numero id da partícula
    pos_vel = pd.concat([n,pos,vel],axis=1)
    pos_vel = pos_vel.sort_values(by=[direction])
    pos_vel.reset_index(drop=True, inplace=True)

    j = int(0)
    concentration_map = np.zeros(int(subspaces))
    velocity_map = np.zeros(int(subspaces))
    if (direction == 'x'):
        aux = dimx/subspaces
        for i in range(int(subspaces)):
            if (j < len(pos_vel)):
                x = pos_vel.loc[j,'x']
            while ((i+1)*aux > x and j < len(pos_vel)-1):
                concentration_map[i] += 1
                velocity_map[i] +=  np.sqrt(pos_vel.loc[j,'v_x']**2 + pos_vel.loc[j,'v_y']**2)
                j += 1 
                x = pos_vel.loc[j,'x']
        mean_Kenergy =  (1/2)*velocity_map**2/concentration_map
        mean_Kenergy[mean_Kenergy == np.inf] = 0
        mean_Kenergy = np.nan_to_num(mean_Kenergy)
        fig, ax1 = plt.subplots()
        ax1.plot(np.linspace(0,dimx,subspaces),mean_Kenergy,'b')
        ax1.set_xlabel("x - direction")
        ax1.set_ylabel("Mean kinetic energy", color='b')

        ax2 = ax1.twinx()
        ax2.plot(np.linspace(0,dimx,subspaces),concentration_map, 'r')
        ax2.set_ylabel("Number of density",color='r')
        #plt.axis([0,dimx,np.min(mean_Kenergy),np.max(mean_Kenergy)])
        fig.tight_layout()
        plt.show(block=True)
    step = input('Choose step or enter -1 to stop:\n')

#def init():
#    concentration_map = [np.nan]*int(subspaces)
#    mean_Kenergy = [np.nan]*int(subspaces)
#    ax1.plot(np.linspace(0,dimx,subspaces),mean_Kenergy,'b')
#    ax1.set_xlabel("x - direction")
#    ax1.set_ylabel("Mean kinetic energy", color='b')
#
#    ax2.plot(np.linspace(0,dimx,subspaces),concentration_map, 'r')
#    ax2.set_ylabel("Number of density",color='r')
#    #plt.axis([0,dimx,np.min(mean_Kenergy),np.max(mean_Kenergy)])
#    fig.tight_layout()    
#
#
#def animate(i):
#    pos = pd.read_csv(res_dir+"position.csv."+step, header=None, names = ["x","y"])
#    vel = pd.read_csv(res_dir+"velocity.csv."+step, header=None, names = ["v_x", "v_y"])
#    n = [x for x in range(len(pos))]
#    n = pd.DataFrame(n, columns=["n"]) # numero id da partícula
#    pos_vel = pd.concat([n,pos,vel],axis=1)
#    pos_vel = pos_vel.sort_values(by=[direction])
#    pos_vel.reset_index(drop=True, inplace=True)
#
#    j = int(0)
#    concentration_map = np.zeros(int(subspaces))
#    velocity_map = np.zeros(int(subspaces))
#    if (direction == 'x'):
#        aux = dimx/subspaces
#        for i in range(int(subspaces)):
#            if (j < len(pos_vel)):
#                x = pos_vel.loc[j,'x']
#            while ((i+1)*aux > x and j < len(pos_vel)-1):
#                concentration_map[i] += 1
#                velocity_map[i] +=  np.sqrt(pos_vel.loc[j,'v_x']**2 + pos_vel.loc[j,'v_y']**2)
#                j += 1 
#                x = pos_vel.loc[j,'x']
#        mean_Kenergy =  (1/2)*velocity_map**2/concentration_map
#        mean_Kenergy[mean_Kenergy == np.inf] = 0
#        mean_Kenergy = np.nan_to_num(mean_Kenergy)
#        ax1.plot(np.linspace(0,dimx,subspaces),mean_Kenergy,'b')
#        ax1.set_xlabel("x - direction")
#        ax1.set_ylabel("Mean kinetic energy", color='b')
#
#        ax2.plot(np.linspace(0,dimx,subspaces),concentration_map, 'r')
#        ax2.set_ylabel("Number of density",color='r')
#        fig.tight_layout()
#        #plt.axis([0,dimx,np.min(mean_Kenergy),np.max(mean_Kenergy)])
#        return ax1 
#        
#        
#
#animate = input('Make animation?(y/n)  ')
#if animate == 'y':
#    fig, ax1 = plt.subplots()
#    ax2 = ax1.twinx()
#    ani = animation.FuncAnimation(fig, animate,frames=n_files, interval=2, save_count=n_files)
#    # ani = animation.FuncAnimation(fig, animate,frames=n_files, init_func=init, interval=2, blit=True, save_count=n_files)
#
#plt.show()