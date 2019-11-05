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
from statsmodels.graphics.tsaplots import plot_acf      

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
    acorr = result[n//2:] / (x.var() * np.arange(n, n//2, -1))
    lag = np.abs(acorr).argmax() + 1
    r = acorr[lag-1]        
#    if np.abs(r) > 0.5:
#      print('Appears to be autocorrelated with r = {}, lag = {}'. format(r, lag))
#    else: 
#      print('Appears to be not autocorrelated')
    return  acorr

# sort the particles into regions
def density_map(x, dimx ,div,ini,fim):
    # input: x: array (1D) of positions in one diraction
    #        dimx: dimension of the region in the diraction of the divison
    #        div: number of divisions
    #        initial position: ini
    #        final position: fin
    # outut: y: list of arrays that contains the indices of the particles that are in a region that is the position in the array
    div = int(dimx/div)
    xp = x // div # This gives the location of each particle
    dmap = [[] for _ in range(div)]
    for i in range(ini,fim):
        dmap[int(xp[i])].append(i)
    return dmap

def comp_mic(x0,x1,mic,dim):
    d = np.array(dim)
    x1 = x1-x0
    mic = mic + (x1 < -d/2)*1
    mic = mic - (x1 > d/2)*1
    return mic 

def cumsum_sq(x):
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = np.sum(np.square([x[0:i]]))
    return y

# Find the directory where the files are
dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
print(dirname)
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

# Read the configuration file
config = configparser.ConfigParser()
config.read(res_dir +  'result/settings.txt')
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])
n_files = int(config['out_files']['out_files'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
t_fim = float(config['global']['t_fim'].split()[0])
dt = float(config['global']['dt'].split()[0])
quant = []
m = []
sigma = []
rs = []
for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
    m.append(float(config['par_'+str(i)]['m'].split()[0]))
    sigma.append(float(config['par_'+str(i)]['sigma'].split()[0])) 
    rs.append(float(config['par_'+str(i)]['rs'].split()[0])) 



# Open the zipped files
zip_rfup = ZipFile(res_dir+'result/rFuP.zip','r')
zip_positions = ZipFile(res_dir+'result/positions.zip','r')
zip_velocities = ZipFile(res_dir+'result/velocities.zip','r')

len_list_files =  len(zip_positions.namelist()) # number of files (steps)

div = 16 #input("Enter in how many regions the region of calculus will be divided in the x diraction [1]: ")
if div == '':
    div = 1
else:
    div = int(div)

Vol = (dimx/div)*dimy
nesb = 100 #int(input("Enter the number of assembles in which this simulation will be divided:\n"))

periodic = True #('y' == input("Were the boundaries periodic? [y/n]\n"))
if periodic:
    len_list_files =  len(zip_rfup.namelist())
mic = 0
ini = 0
fim = quant[0]
if ntype > 0:
    vpg = np.zeros((len_list_files,quant[1],2)) #velocidade das particulas grandes


# len_list_files = 2000
# nesb = 100

esb_len = int(len_list_files/nesb) # length of each ensamble
KE  = np.zeros((len_list_files,div))
densidade = np.zeros((len_list_files,div))
tauk  = np.zeros((len_list_files,div))
tauxyp  = np.zeros((len_list_files,div))
msq  = np.zeros((len_list_files,div))
tauyxp  = np.zeros((len_list_files,div))



print("Reading files...")
if pbar:
    bar = progressbar.ProgressBar(max_value=(len_list_files))
for step in range(1,len_list_files):
    esb = step//esb_len
    if periodic:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["micx","micy","RxFy","RyFx","u","K"])
    else:
        rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","K"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])
    dmap = density_map(pos['x'],dimx,div,ini,fim)

    if ntype > 0:
        vpg[step,:,:] = vel.loc[quant[0]:quant[1]+quant[0]-1,['v_x','v_y']]
    # Depende da região, terá um for

    if step%esb_len == 1: # se step == inicio do ensamble
        if periodic and step > 1:
            r0 = pos[['x','y']] + comp_mic(pos0,pos[['x','y']],mic,[dimx,dimy])*np.array([dimx, dimy]) # um r0 por ensamble
            # r0 = pos[['x','y']] + rfup[["micx", "micy"]].to_numpy()*np.array([dimx, dimy]) # um r0 por ensamble
            lp = dmap # salva as partículas numa daterminada região
        else:
            r0 = pos[['x','y']] # um r0 por ensamble
            lp = dmap # salva as partículas numa daterminada região
    if periodic and step > 1: 
        mic = (comp_mic(pos0, pos[['x','y']],mic,[dimx,dimy])).to_numpy()
        # mic = rfup[["micx", "micy"]].to_numpy()*np.array([dimx, dimy])
    elif step == 1:
        pos0 = pos[['x','y']].to_numpy()
        mic = pos[['x','y']].to_numpy()*0
    
    
    # Separar aqui por região

    for i in range(div):
        #mesmas particulas
        msq[step,i] = np.mean(np.sum(np.square(pos.loc[lp[i],['x','y']] + mic[lp[i],:]*dimy - r0.loc[lp[i],['x','y']]),axis=1))
        densidade[step,i] = len(dmap[i])
        # diferentes particulas
        KE[step,i] = np.mean(rfup.loc[dmap[i],'K'])
        tauk[step,i] = np.sum(vel.loc[dmap[i],'v_x'] * vel.loc[dmap[i],'v_y'] * m[0])/Vol
        tauxyp[step,i] = np.sum(rfup.loc[dmap[i],'RxFy'])/Vol
        tauyxp[step,i] = np.sum(rfup.loc[dmap[i],'RyFx'])/Vol

    pos0 =  pos[['x','y']] #posição anterior

    if pbar: #atualiza a barra de progresso 
        bar.update(step)

kT = KE.mean(axis=0)
plt.figure()
plt.plot(kT)
plt.title("Temperature")

plt.figure()
plt.plot(densidade.mean(axis=0))
plt.title('densidade')

# Calculamos t * viscosidade*2*kT/Vol. kT = KE


etat_kk = np.zeros((esb_len,1))
etat_kp1 = np.zeros((esb_len,1))  
etat_pp1 = np.zeros((esb_len,1))  
etat_kp2 = np.zeros((esb_len,1))  
etat_pp2 = np.zeros((esb_len,1))  
etat_fit1 = np.zeros((esb_len,1))  
etat_fit2 = np.zeros((esb_len,1))  
msq_fit = np.zeros((esb_len,1))  
msq2 = np.zeros((esb_len,1))  

#############  OPÇÃO 1 ########################

if (div > 1):
    overdiv = np.zeros((div,9))

overdiv[:,0] = np.linspace(0,dimx,div)
labels_overdiv = ['x','difu','eta_kk','ETA_XY','eta_xy_kp','eta_xy_pp','ETA_YX','eta_yx_kp','eta_yx_pp','kT','densidade','kn']
for i in range(div):
    etat_kk = 0*etat_kk
    etat_kp1 = 0*etat_kp1
    etat_pp1 = 0*etat_pp1
    etat_kp2 = 0*etat_kp2
    etat_pp2 = 0*etat_pp2
    msq2 = 0*msq2
    for ii in range(1,esb_len): # ao longo de cada ensemble  
        for j in range(nesb): # ao longo dos ensambles
            start = j*esb_len
            fin =  j*esb_len + ii
            etat_kk[ii,0]  += np.trapz(tauk[start:fin,i])**2
            etat_kp1[ii,0] += np.trapz(tauk[start:fin,i])*np.trapz(tauxyp[start:fin,i])
            etat_pp1[ii,0] += np.trapz(tauxyp[start:fin,i])**2 
            etat_kp2[ii,0] += np.trapz(tauk[start:fin,i])*np.trapz(tauyxp[start:fin,i])
            etat_pp2[ii,0] += np.trapz(tauyxp[start:fin,i])**2 

    # Viscosidade * t

    etat_kk  = (etat_kk  * Vol / (2*kT[i]))/nesb 
    etat_kp1 = (etat_kp1 * Vol / (kT[i])  )/nesb 
    etat_pp1 = (etat_pp1 * Vol / (2*kT[i]))/nesb
    etat_kp2 = (etat_kp2 * Vol / (kT[i])  )/nesb 
    etat_pp2 = (etat_pp2 * Vol / (2*kT[i]))/nesb

    # Para descobrir  fazemos um curve fit linear no último terço de pontos
    t = np.linspace(0,esb_len-1, esb_len)*dt
    eta_kk  =  np.polyfit(t[int(esb_len/2):esb_len],etat_kk[int(esb_len/2):esb_len,0],1)
    eta_kp1 = np.polyfit(t[int(esb_len/2):esb_len],etat_kp1[int(esb_len/2):esb_len,0],1)
    eta_pp1 = np.polyfit(t[int(esb_len/2):esb_len],etat_pp1[int(esb_len/2):esb_len,0],1)
    eta_kp2 = np.polyfit(t[int(esb_len/2):esb_len],etat_kp2[int(esb_len/2):esb_len,0],1)
    eta_pp2 = np.polyfit(t[int(esb_len/2):esb_len],etat_pp2[int(esb_len/2):esb_len,0],1)

    # # Calcula a Difusividade
    for j in range(nesb):
        start = j*esb_len
        fin = (j+1)*esb_len
        msq2[:,0] += msq[start:fin,i]
    msq2 = msq2/nesb

    D = np.polyfit(np.arange(0, len(msq2[int(esb_len/3):esb_len,0]),1)*dt,msq2[int(esb_len/3):esb_len],1)
    overdiv[i,1] = D[0]/4

    print("Region: {}\n".format(i))
    fig, axs = plt.subplots(1,3)

    # Plota a difusividade
    msq_fit[:,0] = np.arange(0, len(msq2),1)*dt*D[0] 
    axs[0].scatter(np.arange(0, len(msq2),1) ,msq2)
    axs[0].plot(np.arange(0, len(msq2),1) ,msq_fit,'k')
    axs[0].set_title('Mean square displacement')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    string = "D = {}".format(D[0]/4)
    print(string + "\n")
    axs[0].text(0.95,0.05,string, transform=axs[0].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)

    # Plota as viscosidades usando taupxy

    #plt.figure()
    axs[1].scatter(t,etat_kk,c='g',marker='.',label='t*eta_kk')
    axs[1].scatter(t,etat_kp1,c='b',marker='.',label='t*eta_kp_xy')
    axs[1].scatter(t,etat_pp1,c='r',marker='.',label='t*eta_pp_xy')

    axs[1].plot(t, t*eta_kk[0] + etat_kk[1],'g')
    axs[1].plot(t, t*eta_kp1[0] + etat_kp1[1],'b')
    axs[1].plot(t, t*eta_pp1[0] + etat_pp1[1],'r')

    etat_fit1[:,0] = t*(eta_kk[0] + eta_kp1[0] + eta_pp1[0] ) + etat_kk[1] + etat_kp1[1] + etat_pp1[1]

    axs[1].plot(t,etat_fit1[:,0],'k',linewidth=2, label='t*eta')
    axs[1].legend()
    axs[1].set_title("Viscosity using tau_yx")


    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    string = "etaxy = {:.3}".format(eta_kk[0]+eta_kp1[0]+eta_pp1[0])

    overdiv[i,2] = eta_kk[0]
    overdiv[i,3] = eta_kk[0]+eta_kp1[0]+eta_pp1[0]
    overdiv[i,4] = eta_kp1[0]
    overdiv[i,5] = eta_pp1[0]
    axs[1].text(0.05,0.65,string, transform=axs[1].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)

    # Plota as viscosidades ustando taupyx
    #plt.figure()
    axs[2].scatter(t,etat_kk,c='g',marker='.',label='t*eta_kk')
    axs[2].scatter(t,etat_kp2,c='b',marker='.',label='t*eta_kp_yx')
    axs[2].scatter(t,etat_pp2,c='r',marker='.',label='t*eta_pp_yx')

    axs[2].plot(t, t*eta_kk[0] + etat_kk[1],'g')
    axs[2].plot(t, t*eta_kp2[0] + etat_kp2[1],'b')
    axs[2].plot(t, t*eta_pp2[0] + etat_pp2[1],'r')

    etat_fit2[:,0] = t*(eta_kk[0] + eta_kp2[0] + eta_pp2[0] ) + etat_kk[1] + etat_kp2[1] + etat_pp2[1]
    axs[2].plot(t, t*(eta_kk[0] + eta_kp2[0] + eta_pp2[0] ) + etat_kk[1] + etat_kp2[1] + etat_pp2[1],'k',linewidth=2, label='t*eta')
    axs[2].legend()
    axs[2].set_title("Viscosity using tau_xy.")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    string = "etayx = {:.3f}".format(eta_kk[0]+eta_kp2[0]+eta_pp2[0])
    overdiv[i,6] = eta_kk[0]+eta_kp2[0]+eta_pp2[0]
    overdiv[i,7] = eta_kp2[0]
    overdiv[i,8] = eta_pp2[0]
    axs[2].text(0.05,0.65,string, transform=axs[2].transAxes, fontsize=10, verticalalignment='bottom', bbox=props)

    print("\nCom tau_xy: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[0],eta_kp1[0],eta_pp1[0], eta_kk[0]+eta_kp1[0]+eta_pp1[0]))
    print("\nCom tau_x]yx: \neta_kk = {}\neta_kp = {}\neta_pp = {}\neta = {}\n".format(eta_kk[0],eta_kp2[0],eta_pp2[0], eta_kk[0]+eta_kp2[0]+eta_pp2[0]))

    fig.suptitle("Region {}".format(i))
    t = t.reshape((len(t),1))
    csv_out = pd.DataFrame(data=np.concatenate((etat_kk,etat_kp1,etat_pp1,etat_fit1,etat_kp2,etat_pp2, etat_fit2, msq2, msq_fit,t),axis=1), columns=["eta_kk","eta_kp_xy","eta_pp_xy","fit_eta_xy","eta_kp_yx","eta_pp_yx","fit_eta_yx","msq","msq_fit","t"])

    csv_out.to_csv(res_dir + "Region_{}.csv".format(i),index=False,sep=' ')
plt.figure()
plt.plot(overdiv[:,3],label='eta_xy')
plt.plot(overdiv[:,6],label='eta_yx')
plt.legend()
plt.title('eta Over the positions')

plt.figure()
plt.plot(overdiv[:,2],label='eta_kk')
plt.plot(overdiv[:,4],label='eta_kp')
plt.plot(overdiv[:,5],label='eta_pp')
plt.legend()
plt.title('eta Over the positions')

plt.figure()
plt.plot(overdiv[:,1],label='D')
plt.title('Dif Over the positions')

densidade = densidade.mean(axis=0)/Vol
densidade = np.reshape(densidade,(div,1))
kT = np.reshape(kT, (div,1))
mfp = 1/(2**.5 * 2*sigma[0] * densidade)
kn = mfp/5 
plt.figure()
plt.plot(kn)
plt.title("kn")


odv = pd.DataFrame(data=np.hstack((overdiv,kT,densidade,kn)),columns=labels_overdiv)

odv.to_csv(res_dir + 'over_regions.csv',index=False,sep=' ')

xpg = pos.loc[quant[0]:quant[1]+quant[0]-1,['x','y']].to_numpy()   

select = []
# seleciona só as partículas que estão 10% longe da parede
for i in range(len(xpg)):
    if abs(xpg[i,0]-dimx/2) < (dimx/2)*0.8:
        select.append(i)


acrx = []
acry = []
Rx = []
Ry = []
R = []
acr = []
a = input('Number of steps for the Autocorrelation ')
if a == '':
    vpg2 = vpg
else:
    a = int(a)
    vpg2 = np.zeros((a,vpg.shape[1],vpg.shape[2]))
    j = 0
    for i in range(len(vpg)):
        b = int(len(vpg)/a)
        if int(i%b) == 0 and j < len(vpg2):
            vpg2[j,:,:] = vpg[i,:,:]
            j += 1

for i in select:
    acr1 = autocorr(vpg2[:,i,0])
    Rx.append(np.trapz(acr1))
    acrx.append(acr1)

    acr2 = autocorr(vpg2[:,i,1])
    Ry.append(np.trapz(acr2))
    acry.append(acr2)''

    acr12 = autocorr(np.sqrt(vpg2[:,i,1]**2 + vpg2[:,i,0]**2 ))
    R.append(np.trapz(acr12))
    acr.append(acr12)

Rx = np.array(Rx)
Ry = np.array(Ry)
a = pos.loc[[i for i in range(quant[0],quant[0]+quant[1])],'x'].to_numpy()
plt.figure()
plt.scatter(a[select],Rx)
acrx = np.array(acrx)
acry = np.array(acry)
acr = np.array(acr)
plt.figure()

meanacr1 = np.mean(acrx,axis=0)
conf = np.arange(len(meanacr1))
conf = 1.96/(np.sqrt(-conf + len(vpg2)))
conf_arima = 1.96*(np.sqrt( (1 + 2*cumsum_sq(meanacr1))/(-conf + len(vpg2)) ))
plt.plot(meanacr1, label="vx")
meanacr2 = np.mean(acry,axis=0)
plt.plot(meanacr2, label="vy")
meanacr3 = np.mean(acr,axis=0)
plt.plot(meanacr3, label="v")
plt.plot(conf,'.r')
plt.plot(-conf,'.r')
plt.plot(conf_arima,'-r')
plt.plot(-conf_arima,'-r')

plt.legend()

meanacr1 = meanacr1.reshape((len(meanacr1),1))
meanacr2 = meanacr2.reshape((len(meanacr2),1))
meanacr3 = meanacr3.reshape((len(meanacr3),1))
conf = conf.reshape((len(conf),1))
conf_arima = conf_arima.reshape((len(conf_arima),1))

df2 = pd.DataFrame(data=np.hstack([meanacr1, meanacr2, meanacr3,conf,-conf,conf_arima,-conf_arima]),columns=["axrcx","acrcy","acr","conf1","conf2","conf_arima","cong_arima2"])

df2.to_csv(res_dir + 'acr.csv',index=False,sep=' ')

plt.show()
