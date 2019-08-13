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

def einstein_relation(A,na):
    itv = int(len(A)/na)
    B = np.zeros((itv,1))
    for j in range(itv):
        for i in range(na):
            B[j] += (A[i*itv+j] - A[i*itv])**2
    B = B/na 
    return B
    
    

dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

zip_rfup = ZipFile(res_dir+'/rFuP.zip','r')
zip_positions = ZipFile(res_dir+'/positions.zip','r')
zip_velocities = ZipFile(res_dir+'/velocities.zip','r')

config = configparser.ConfigParser()
config.read(res_dir + 'settings.txt')
N = int(config['global']['N'].split()[0])
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])
n_files = int(config['out_files']['out_files'].split()[0])

ntype = int(config['global']['Ntype'].split()[0])
t_fim = float(config['global']['t_fim'].split()[0])
dt = float(config['global']['dt'].split()[0])
F = np.array([0,0])
quant = []
sigma = []
epsilon = []
rs = [] # raio sólido
mass = []
tipo = [0]*N

dt = t_fim/n_files

a = input("Enter the subdomains mesh dimensions.\n")
a = a.split()
mesh = np.array([int(a[0]), int(a[1])])
a = input("Enter a location. xmin xmax ymin ymax\n")
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

nsteps = n_files # int(input('Enter the number of steps:\n'))

density_map = np.zeros((mesh[0],mesh[1], nsteps+1))
r = np.zeros((mesh[0],mesh[1],nsteps+1))
dQm = np.zeros(nsteps)
eta = np.zeros(nsteps)


step = 0 
n1,n2 = 0,0
stdscr = 's' #para barra de progresso
sample_list = []
Qxy = np.zeros((nsteps+1,1))
Qyx = np.zeros((nsteps+1,1))
Npart = np.zeros((nsteps+1,1))
KE = np.zeros((nsteps+1,1)) # kT = KE
#eta = np.zeros((nsteps,1))

if pbar:
    bar = progressbar.ProgressBar(max_value=nsteps)
while step <= nsteps: 
    r0 = []
    p0 = []

    particle_map = [[[] for _ in range(mesh[1])] for _ in range(mesh[0])]
    rfup = pd.read_csv(zip_rfup.open('rF_u_P.csv.'+str(step)), header=None, names = ["RxFy","RyFx","u","px","py"])
    pos = pd.read_csv(zip_positions.open("position.csv."+str(step)), header=None, names = ["x","y"])
    vel = pd.read_csv(zip_velocities.open("velocity.csv."+str(step)), header=None, names = ["v_x", "v_y"])

    n = [x for x in range(len(pos))]
    n = pd.DataFrame(n, columns=["n"]) # numero id da partícula
    pos_vel = pd.concat([n,pos,vel,tipo],axis=1)
   
    for nn in range(len(pos_vel)):
        xp = int(pos_vel.loc[nn,'x']//hx)
        yp = int(pos_vel.loc[nn,'y']//hy)
        if xp == mesh[0]:
            xp = xp - 1
        if yp == mesh[1]:
            yp = yp - 1
        particle_map[xp][yp].append( pos_vel.loc[nn,'n'] )
        density_map[xp,yp,step] += 1


    for i in range(region[0],region[1]):
        for j in range(region[2],region[3]):
            nparcl = len(particle_map[i][j])
            for nn in range(nparcl):
                Npart[step] += 1 #numero de partículas na região neste passo 
                n1 = particle_map[i][j][nn]
                r0.append([pos_vel.loc[n1,'x'], pos_vel.loc[n1,'y']])
                m = mass[pos_vel.loc[n1,'tipo']]
                p0.append([pos_vel.loc[n1,'px'], pos_vel.loc[n1,'py']])
                KE[step] += np.sqrt(pos_vel.loc[n1,'v_x']**2 + pos_vel.loc[n1,'v_y']**2)*m/nparcl
   
    r_0 = np.array(r0)
    p_0 = np.array(p0)
    Qxy[step] += np.sum(r_0[:,0]*p_0[:,1], axis=0)/Vol 
    Qyx[step] += np.sum(r_0[:,1]*p_0[:,0], axis=0)/Vol
    if pbar:
        bar.update(step)
    step += 1

KE2 = np.sum(KE[1:])/len(KE[1:])
na = int(input('Enter de number of assembles (na) to average.\nThe number of steps for correlation will be nsteps/na: '))

etaxyt = einstein_relation(Qxy,na)*Vol/(2*KE2)
etayxt = einstein_relation(Qyx,na)*Vol/(2*KE2)

aux = len(etaxyt)
t = np.linspace(0,t_fim/na,aux)
plt.figure(1)
plt.plot(t,etaxyt,'.r', label='etaxy*t')
m,b = np.polyfit(t, etaxyt, 1)
y = m*t + b 
plt.plot(t,y,'-r', label='etaxy*t fit')
plt.plot(t,etayxt,'.b', label='etayx')
m,b = np.polyfit(t[int(
    /2):], etayxt[int(aux/2):], 1)
y = m*t + b 
plt.plot(t,y,'-b', label='etaxy*t fit')
print("eta = {}".format(m[0]))
plt.legend()

# h = t_fim/nsteps 
# etaxy = diff1(etaxy2t,h)/(2*KE2)
# etayx = diff1(etayx2t,h)/(2*KE2)
# plt.figure(2)
# plt.plot(etaxy[1:len(etaxy)], label='etaxy*kT')
# plt.plot(etayx[1:len(etaxy)], label='etayx*kT')
# plt.plot((etayx[1:len(etaxy)]+etaxy[1:len(etaxy)])/2,'--k', label='eta*kT mean')
# plt.legend()

plt.show()