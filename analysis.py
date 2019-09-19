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
    div = int(dimx/div)
    xp = x // div # This gives the location of each particle
    dmap = [[] for _ in range(div)]
    for i in range(len(x)):
        dmap[int(xp[i])].append(i)
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

KE  = np.zeros((len_list_files,div))
tauk  = np.zeros((len_list_files,div))
tauxyp  = np.zeros((len_list_files,div))
msq  = np.zeros((len_list_files,div))
tauyxp  = np.zeros((len_list_files,div))

nesb = int(input("Enter the number of assembles in which this simulation will be divided:\n"))
periodic = ('y' == input("Were the boundaries periodic? [y/n]\n"))
esb_len = int(len_list_files/nesb) # length of each ensamble
mic = 0
print("Reading files...")
if pbar:
    bar = progressbar.ProgressBar(max_value=(len_list_files))
for step in range(1,len_list_files):
    esb = step//esb_len
    if periodic:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["micx","micy","RxFy","RyFx","u","m","K"])
    else:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","m","K"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])

    dmap = density_map(pos['x'],dimx,div)

    # Depende da região, terá um for
    if step%esb_len == 1: # se step == inicio do ensamble
        r0 = pos[['x','y']] # um r0 por ensamble
        lp = dmap # salva as partículas numa daterminada região
    # if periodic: 
        # mic = rfup[["micx", "micy"]].to_numpy()*np.array([dimx, dimy])
    
    # Separar aqui por região

    for i in range(div):
        #mesmas particulas
        # msq[step,i] = np.mean(np.sum(np.square(pos.loc[lp[i],['x','y']] + mic[lp[i],:] - r0.loc[lp[i]]),axis=1))
        # diferentes particulas
        KE[step,i] = np.mean(rfup.loc[dmap[i],'K'])
        tauk[step,i] = np.sum(vel.loc[dmap[i],'v_x'] * vel.loc[dmap[i],'v_y'] * rfup.loc[dmap[i],'m'])/Vol
        tauxyp[step,i] = np.sum(rfup.loc[dmap[i],'RxFy'])/Vol
        tauyxp[step,i] = np.sum(rfup.loc[dmap[i],'RyFx'])/Vol

    if pbar: #atualiza a barra de progresso 
        bar.update(step)

kT = KE.mean(axis=0)

# Calculamos t * viscosidade*2*kT/Vol. kT = KE


etat_kk = np.zeros(esb_len)
etat_kp1 = np.zeros(esb_len)
etat_pp1 = np.zeros(esb_len)
etat_kp2 = np.zeros(esb_len) 
etat_pp2 = np.zeros(esb_len)

if (div > 1):
    overdiv = np.zeros((div,2))

for i in range(div):
    etat_kk = 0*etat_kk
    etat_kp1 = 0*etat_kp1
    etat_pp1 = 0*etat_pp1
    etat_kp2 = 0*etat_kp2
    etat_pp2 = 0*etat_pp2
    for ii in range(1,esb_len): # ao longo de cada ensemble  
        for j in range(nesb): # ao longo dos ensambles
            start = j*esb_len
            fin =  j*esb_len + ii
            etat_kk[ii]  += np.trapz(tauk[start:fin,i])**2
            etat_kp1[ii] += np.trapz(tauk[start:fin,i])*np.trapz(tauxyp[start:fin,i])
            etat_pp1[ii] += np.trapz(tauxyp[start:fin,i])**2 
            etat_kp2[ii] += np.trapz(tauk[start:fin,i])*np.trapz(tauyxp[start:fin,i])
            etat_pp2[ii] += np.trapz(tauyxp[start:fin,i])**2 

    # Viscosidade * t

    etat_kk  = (etat_kk  * Vol / (2*kT[i]))/nesb 
    etat_kp1 = (etat_kp1 * Vol / (kT[i])  )/nesb 
    etat_pp1 = (etat_pp1 * Vol / (2*kT[i]))/nesb
    etat_kp2 = (etat_kp2 * Vol / (kT[i])  )/nesb 
    etat_pp2 = (etat_pp2 * Vol / (2*kT[i]))/nesb

    # Para descobrir  fazemos um curve fit linear no último terço de pontos
    t = np.linspace(0,esb_len-1, esb_len)
    eta_kk = np.polyfit(t[int(esb_len/2):esb_len],etat_kk[int(esb_len/2):esb_len],1)
    eta_kp1 = np.polyfit(t[int(esb_len/2):esb_len],etat_kp1[int(esb_len/2):esb_len],1)
    eta_pp1 = np.polyfit(t[int(esb_len/2):esb_len],etat_pp1[int(esb_len/2):esb_len],1)
    eta_kp2 = np.polyfit(t[int(esb_len/2):esb_len],etat_kp2[int(esb_len/2):esb_len],1)
    eta_pp2 = np.polyfit(t[int(esb_len/2):esb_len],etat_pp2[int(esb_len/2):esb_len],1)

    # # Calcula a Difusividade
    # D = np.polyfit(t,msq,1)

    print("Region: {}\n".format(i))
    fig, axs = plt.subplots(1,3)
    # axs[1,0].figure()

    # # Plota a difusividade
    # axs[1,0].scatter(t,msq[:,i])
    # axs[1,0].plot(t,t*D[0] + D[1],'k')
    # axs[1,0].title('Mean square displacement')
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # string = "D = {}".format(D[0]/4)
    # print(string + "\n")
    # axs[1,0].text(0.95,0.05,string, transform=axs[1,0].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)


    # Plota as viscosidades usando taupxy

    #plt.figure()
    axs[1].scatter(t,etat_kk,c='g',marker='.',label='t*eta_kk')
    axs[1].scatter(t,etat_kp1,c='b',marker='.',label='t*eta_kp_xy')
    axs[1].scatter(t,etat_pp1,c='r',marker='.',label='t*eta_pp_xy')

    axs[1].plot(t, t*eta_kk[0] + etat_kk[1],'g')
    axs[1].plot(t, t*eta_kp1[0] + etat_kp1[1],'b')
    axs[1].plot(t, t*eta_pp1[0] + etat_pp1[1],'r')

    axs[1].plot(t, t*(eta_kk[0] + eta_kp1[0] + eta_pp1[0] ) + etat_kk[1] + etat_kp1[1] + etat_pp1[1],'k',linewidth=2, label='t*eta')
    axs[1].legend()
    axs[1].set_title("Viscosity using tau_yx")

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    string = "etaxy = {:.3}".format(eta_kk[0]+eta_kp1[0]+eta_pp1[0])
    overdiv[i,0] = eta_kk[0]+eta_kp1[0]+eta_pp1[0]
    axs[1].text(0.05,0.65,string, transform=axs[1].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)

    # Plota as viscosidades ustando taupyx
    #plt.figure()
    axs[2].scatter(t,etat_kk,c='g',marker='.',label='t*eta_kk')
    axs[2].scatter(t,etat_kp2,c='b',marker='.',label='t*eta_kp_yx')
    axs[2].scatter(t,etat_pp2,c='r',marker='.',label='t*eta_pp_yx')

    axs[2].plot(t, t*eta_kk[0] + etat_kk[1],'g')
    axs[2].plot(t, t*eta_kp2[0] + etat_kp2[1],'b')
    axs[2].plot(t, t*eta_pp2[0] + etat_pp2[1],'r')

    axs[2].plot(t, t*(eta_kk[0] + eta_kp2[0] + eta_pp2[0] ) + etat_kk[1] + etat_kp2[1] + etat_pp2[1],'k',linewidth=2, label='t*eta')
    axs[2].legend()
    axs[2].set_title("Viscosity using tau_xy.")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    string = "etayx = {:.3f}".format(eta_kk[0]+eta_kp2[0]+eta_pp2[0])
    overdiv[i,1] = eta_kk[0]+eta_kp2[0]+eta_pp2[0]
    axs[2].text(0.05,0.65,string, transform=axs[2].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)

    print("\nCom tau_xy: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[0],eta_kp1[0],eta_pp1[0], eta_kk[0]+eta_kp1[0]+eta_pp1[0]))
    print("\nCom tau_x]yx: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[0],eta_kp2[0],eta_pp2[0], eta_kk[0]+eta_kp2[0]+eta_pp2[0]))

    fig.suptitle("Region {}".format(i))

plt.figure()
plt.plot(overdiv[:,0],label='eta_xy')
plt.plot(overdiv[:,1],label='eta_yx')
plt.title('Over the positions')

plt.show()