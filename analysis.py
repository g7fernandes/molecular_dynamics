''' This program analyses viscosity, kn, pressure and other properties of a simulation as a whole
The input must be zipped files of position, velocity, and potential and force times tranversal displacement'''

import numpy as np 
import pandas as pd 
import os
import configparser
import matplotlib.pyplot as plt
from glob import glob
import time
from zipfile import ZipFile


try:
    import progressbar
    pbar = True 
except:
    print('progressbar not available. Try one:')
    print('conda install -c conda-forge progressbar2')
    print('conda install -c conda-forge/label/gcc7 progressbar2')
    pbar = False


# Find the directory where the files are
dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

# Read the configuration file
config = configparser.ConfigParser()
config.read(res_dir + '/settings.txt')
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])
Vol = dimx*dimy

# Open the zipped files
zip_rfup = ZipFile(res_dir+'/rFuP.zip','r')
zip_positions = ZipFile(res_dir+'/positions.zip','r')
zip_velocities = ZipFile(res_dir+'/velocities.zip','r')

len_list_files =  len(zip_positions.namelist()) # number of files (steps)

KE = np.zeros(len_list_files)
tauk = np.zeros(len_list_files)
tauxyp = np.zeros(len_list_files)
tauyxp = np.zeros(len_list_files)
msq = np.zeros(len_list_files)

nesb = int(input("Enter the number of assembles in which this simulation will be divided:\n"))
periodic = ('y' == input("Where the boundaries periodic? [y/n]\n"))
for step in range(len_list_files):
    if periodic:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["micx","micy","RxFy","RyFx","u","m","K"])
    else:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","m","K"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])

    if step == 0:
        r0 = pos[['x','y']]

    KE[step] = np.sum(rfup('K'))
    tauk[step] = np.sum(vel['v_x'] * vel['v_y'] * rfup ['m'])/Vol
    tauxyp[step] = np.sum('RxFy')/Vol
    tauyxp[step] = np.sum('RyFx')/Vol
    
    if periodic: 
        mic = rfup[["micx", "micy"]]*np.array([dimx, dimy])
        msq[step] = np.sum(np.sum(np.square(pos[['x','y']] + mic - r0),axis=1))
    
# Difusividade 2D
D = msq/(4*(len(r0)-1))

esb_len = int(len_list_files)/nesb # length of each assemble
for i in range(1,esb_len): # ao longo de cada ensemble  
    for j in range(nesb): # ao longo dos ensambles
        start = j*esb_len
        fin =  j*esb_len + i
        eta_kk = np.trapz(tauk[start:fin])**2
        
        eta_kp1 = np.trapz(tauk[start:fin])*np.trapz(tauxyp[start:fin])
        eta_pp1 = np.trapz(tauxyp[start:fin])**2 

        eta_kp2 = np.trapz(tauk[start:fin])*np.trapz(tauyxp[start:fin])
        eta_pp2 = np.trapz(tauyxp[start:fin])**2 
    
