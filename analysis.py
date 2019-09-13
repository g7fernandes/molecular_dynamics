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

# Autocorrelation function
def autocorr(x):
    n = x.size
    norm = (x - np.mean(x))
    result = np.correlate(norm, norm, mode='same')
    acorr = result[n//2 + 1:] / (x.var() * np.arange(n-1, n//2, -1))
    lag = np.abs(acorr).argmax() + 1
    r = acorr[lag-1]        
    if np.abs(r) > 0.5:
      print('Appears to be autocorrelated with r = {}, lag = {}'. format(r, lag))
    else: 
      print('Appears to be not autocorrelated')
    return r, lag

# sort the particles into regions
def density_map(x, dimx ,div):
    # input: x: array (1D) of positions in one diraction
    #        dimx: dimension of the region in the diraction of the divison
    #        div: number of divisions
    # outut: y: list of arrays that contains the indices of the particles that are in a region that is the position in the array
    div = dimx/div
    xp = x // div # This gives the location of each particle
    dmap = [[] for _ in range(div)]
    for i in range(len(x)):
        dmap[xp].append(i)
    return dmap

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

div = input("Enter in how many regions the region of calculus will be divided in the x diraction [1]: ")
if div == '':
    div = 1
else:
    div = int(div)

KE,tauk,tauxyp,msq,tauyxp  = np.zeros((len_list_files,div))

nesb = int(input("Enter the number of assembles in which this simulation will be divided:\n"))
periodic = ('y' == input("Where the boundaries periodic? [y/n]\n"))
mic = 0
for step in range(len_list_files):
    if periodic:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["micx","micy","RxFy","RyFx","u","m","K"])
    else:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","m","K"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])

    dmap = density_map(pos['x'],dimx,div)

    # Depende da região, terá um for
    if step == 0: # se step == inicio do ensamble
        r0 = pos[['x','y']] # um r0 por ensamble
        # lp = dmap # salva as partículas numa daterminada região
    if periodic: 
        mic = rfup[["micx", "micy"]]*np.array([dimx, dimy])
    
    # Separar aqui por região
    msq[step] = np.sum(np.sum(np.square(pos[['x','y']] + mic - r0),axis=1))

    KE[step] = np.sum(rfup('K'))
    tauk[step] = np.sum(vel['v_x'] * vel['v_y'] * rfup ['m'])/Vol
    tauxyp[step] = np.sum('RxFy')/Vol
    tauyxp[step] = np.sum('RyFx')/Vol
    

    
# Difusividade 2D
D = msq/(4*(len(r0)-1))

kT = np.mean(KE)

# Calculamos t * viscosidade*2*kT/Vol. kT = KE

esb_len = int(len_list_files)/nesb # length of each ensamble
etat_kk, etat_kp1, etat_pp1, etat_kp2, etat_pp2 = np.zeros(esb_len)


for i in range(1,esb_len): # ao longo de cada ensemble  
    for j in range(nesb): # ao longo dos ensambles
        start = j*esb_len
        fin =  j*esb_len + i
        etat_kk[i]  += np.trapz(tauk[start:fin])**2
        etat_kp1[i] += np.trapz(tauk[start:fin])*np.trapz(tauxyp[start:fin])
        etat_pp1[i] += np.trapz(tauxyp[start:fin])**2 
        etat_kp2[i] += np.trapz(tauk[start:fin])*np.trapz(tauyxp[start:fin])
        etat_pp2[i] += np.trapz(tauyxp[start:fin])**2 

# Viscosidade * t

etat_kk  = (etat_kk  * Vol / (2*kT))/nesb 
etat_kp1 = (etat_kp1 * Vol / (kT)  )/nesb 
etat_pp1 = (etat_pp1 * Vol / (2*kT))/nesb
etat_kp2 = (etat_kp2 * Vol / (kT)  )/nesb 
etat_pp2 = (etat_pp2 * Vol / (2*kT))/nesb

# Para descobir T fazemos um curve fit linear 
t = np.linspace(0,esb_len-1, esb_len)
eta_kk = np.polyfit(t,etat_kk,1)
eta_kp1 = np.polyfit(t,etat_kp1,1)
eta_pp1 = np.polyfit(t,etat_pp1,1)
eta_kp2 = np.polyfit(t,etat_kp2,1)
eta_pp2 = np.polyfit(t,etat_pp2,1)

# Plota as viscosidades usando taupxy

plt.figure()
plt.scatter(t,etat_kk,c='g',label='t*eta_kk')
plt.scatter(t,etat_kp1,c='b',label='t*eta_kp_xy')
plt.scatter(t,etat_pp1,c='r',label='t*eta_pp_xy')

plt.plot(t, t*eta_kk[1] + etat_kk[0],'g')
plt.plot(t, t*eta_kp1[1] + etat_kp1[0],'b')
plt.plot(t, t*eta_pp1[1] + etat_pp1[0],'r')

plt.plot(t, t*(eta_kk[1] + eta_kp1[1] + eta_pp1[1] ) + etat_kk[0] + etat_kp1[0] + etat_pp1[0],'k',linewidth=2, label='t*eta')
plt.legend()
plt.title("Visocisdade usando tau_yx")
# Plota as viscosidades ustando taupyx
plt.figure()
plt.scatter(t,etat_kk,c='g',label='t*eta_kk')
plt.scatter(t,etat_kp2,c='b',label='t*eta_kp_yx')
plt.scatter(t,etat_pp2,c='r',label='t*eta_pp_yx')

plt.plot(t, t*eta_kk[1] + etat_kk[0],'g')
plt.plot(t, t*eta_kp2[1] + etat_kp2[0],'b')
plt.plot(t, t*eta_pp2[1] + etat_pp2[0],'r')

plt.plot(t, t*(eta_kk[1] + eta_kp2[1] + eta_pp2[1] ) + etat_kk[0] + etat_kp2[0] + etat_pp2[0],'k',linewidth=2, label='t*eta')
plt.legend()
plt.title("Visocisdade usando tau_xy")
print("\nCom tau_xy: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[1],eta_kp1[1],eta_pp1[1], eta_kk[1]+eta_kp1[1]+eta_pp1[1]))
print("\nCom tau_x]yx: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[1],eta_kp2[1],eta_pp2[1], eta_kk[1]+eta_kp2[1]+eta_pp2[1]))

plt.show()