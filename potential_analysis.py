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

kB = 1.38064852*10**-23

def integral(F,i,j,dt):
    I = np.trapz(F[i,j,:],dx=dt)                
    return I

def diff1(F,h):
    B = np.zeros((len(F),1))
    for i in range(len(F)):
        if i > 1 and i < len(F) - 2:
            B[i] = ((1/12)*F[i-2] + (-2/3)*F[i-1] + (2/3)*F[i+1] + (-1/12)*F[i+2])/h   
        elif i == 0:
            B[i] = ((-11/6)*F[0] + 3*F[1] + (-3/2)*F[2] + (1/3)*F[3])/h 
        elif i == 1:
            B[i] = (-2*F[0]-3*F[1]+6*F[2]-1*F[3])/(6*h)
        elif i == len(F)-1:
            B[i] = (1*F[i-2]-4*F[i-1]+3*F[i+0])/(2*1.0*h**1)
        elif i == len(F)-2:
            B[i] = (1*F[i-2]-6*F[i-1]+3*F[i+0]+2*F[i+1])/(6*1.0*h**1)
    return B

def greenkubo(Tk,Tp,i,j,na,dt):
    itv = int(Tp.shape[2]/na)
    kk = np.zeros((itv))
    kp = np.zeros((itv))
    pp = np.zeros((itv))
    for a in range(itv):
        for b in range(na):
            kk[a] += np.trapz(Tk[i,j,b*itv:b*itv+a],dx=dt)**2
            kp[a] += np.trapz(Tp[i,j,b*itv:b*itv+a],dx=dt)*np.trapz(Tk[i,j,b*itv:b*itv+a],dx=dt)
            pp[a] += np.trapz(Tp[i,j,b*itv:b*itv+a],dx=dt)**2
    kk = kk/na
    pp = pp/na
    kp = kp/na 
    return kk, kp, pp

def einstein_condutividade(delta_e,na):
    lambda_T = np.zeros(delta_e.shape[0], delta_e.shape[1], na)
    itv = int(delta_e.shape[2]/na)
    for a  in range(itv):
        for b in range(na):
            for i in range(len(delta_e.shape[0])):
                for j in range(len(delta_e.shape[1])):
                    lambda_T[i,j,a] += (delta_e[i,j,b*itv] - delta_e[i,j,b*itv+a])**2
    lambda_T = lambda_T/na 
    return lambda_T

    
dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

#zip_rfup = ZipFile(res_dir+'/rFuP.zip','r')
zip_positions = ZipFile(res_dir+'/positions.zip','r')
zip_velocities = ZipFile(res_dir+'/velocities.zip','r')

len_list_files =  len(zip_positions.namelist()+zip_velocities.namelist())

config = configparser.ConfigParser()
config.read(res_dir + 'settings.txt')
N = int(config['global']['N'].split()[0])
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])
n_files = int(config['out_files']['out_files'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
t_fim = float(config['global']['t_fim'].split()[0])
#dt = float(config['global']['dt'].split()[0])
F = np.array([0,0])
quant = []
sigma = []
epsilon = []
rs = [] # raio solido
mass = []
tipo = [0]*N

dt = t_fim/n_files

a = input("Enter the subdomains mesh dimensions.\n")
a = a.split()
mesh = np.array([int(a[0]), int(a[1])])
a = input("Enter a location (Starts at 0). xmin xmax ymin ymax\n")
if a == '':
    region = [0, mesh[0], 0, mesh[1]]
else:
    region = [int(x) for x in a.split()]

Vol = (dimx/mesh[0])*(dimy/mesh[1]) # volume dos elementos da malha
for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
    rs.append(float(config['par_'+str(i)]['rs'].split()[0]))
    sigma.append(float(config['par_'+str(i)]['sigma'].split()[0]))
    epsilon.append(float(config['par_'+str(i)]['epsilon'].split()[0]))
    mass.append(float(config['par_'+str(i)]['m'].split()[0]))
j,k = 0,0
for i in range(len(quant)):
    for j in range(quant[i]):
        tipo[j+k] = i
        
    k = sum(quant[0:i+1])
tipo = pd.DataFrame(tipo, columns=["tipo"]) # numero id da partícula

hx = dimx/mesh[0]
hy = dimy/mesh[1]

print("Planned output {} files (steps).\n".format(n_files))
nsteps = int(input('Enter the number of steps (final number):\n'))


print("There are {} files.".format(n_files-1))
try:
    step = int(input("Enter the initial step [0]: "))
except:
    step = 0
    
density_map = np.zeros((mesh[0],mesh[1], nsteps+1 - step))
Kmap = np.zeros((mesh[0],mesh[1],nsteps+1 - step))
Vmap = np.zeros((mesh[0],mesh[1],nsteps+1 - step))
a = nsteps - step
stepini = step 
if pbar:
    bar = progressbar.ProgressBar(max_value=(a))
while step <= nsteps: 
    particle_map = [[[] for _ in range(mesh[1])] for _ in range(mesh[0])]
#    rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","px","py"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])
    n = [x for x in range(len(pos))]
    n = pd.DataFrame(n, columns=["n"]) # numero id da partícula
    # print("AA")
    pos_vel = pd.concat([n,pos,vel,tipo],axis=1)
    # print("B")
    for nn in range(quant[0]): # só vai contar as particulas do grupo 1
        xp = int(pos_vel.loc[nn,'x']//hx)
        yp = int(pos_vel.loc[nn,'y']//hy)
        if xp == mesh[0]:
            xp = xp - 1
        if yp == mesh[1]:
            yp = yp - 1
        particle_map[xp][yp].append( pos_vel.loc[nn,'n'] )
        density_map[xp,yp,step-stepini] += 1
    # "print C"
    for i in range(region[0],region[1]):
        for j in range(region[2],region[3]):
            for nn in range(len(particle_map[i][j])):
                n1 = particle_map[i][j][nn]
                m = mass[pos_vel.loc[n1,'tipo']]
                Kmap[i,j,step-stepini] += (pos_vel.loc[n1,'v_x']**2+pos_vel.loc[n1,'v_y']**2)*m/(2)
                # Vmap[i,j,step] += pos_vel.loc[n1,'u']     
    if pbar:
        bar.update(step-stepini)
    step += 1

KE = Kmap/density_map

KE = np.nan_to_num(KE)
KE = np.sum(KE,axis=2)/a
plt.figure()
if mesh[1] == 1:
    T = KE #np.sum(KE,axis=0) # energia cinética ao longo de X
    plt.plot(T,label='Temperatura adimensional ao longo de x')
elif mesh[0] == 1:
    T = KE #np.sum(KE,axis=1) # energia cinética ao longo de Y
    plt.plot(T,label='Temperatura adimensional ao longo de y')
plt.legend()

# plt.figure()

# Vlj = Vmap/a

# if mesh[1] == 1:
#     V = np.sum(Vlj,axis=0) # energia cinética ao longo de X
#     plt.plot(T,label='Potencial adimensional ao longo de x')
# elif mesh[0] == 1:
#     V = np.sum(Vlj,axis=1) # energia cinética ao longo de Y
#     plt.plot(T,label='Potencial adimensional ao longo de y')
# plt.legend()


plt.show()


