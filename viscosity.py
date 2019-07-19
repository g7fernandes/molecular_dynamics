import numpy as np 
import pandas as pd 
import os
import configparser
import matplotlib.pyplot as plt
from glob import glob
import time
try:
    import progressbar
    pbar = True 
except:
    print('progressbar not available. Try one:')
    print('conda install -c conda-forge progressbar2')
    print('conda install -c conda-forge/label/gcc7 progressbar2')
    pbar = False

def integral(F,i,j,dt):
    # sh = np.shape(F)
    I = np.trapz(F[i,j,:],dx=dt)                
    # I = np.zeros((sh[0],sh[1]))
    # xi = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
    #xp = np.linspace(0, dt*sh[2], sh[2])
    #w = [5/9, 8/9, 5/9]
    # for i in range(sh[0]):
    #     for j in range(sh[1]):
            # for s in range(sh[2]-1):
                # I[i,j] += 0.5*(dt)*( \
                #     w[0]*np.interp( dt*xi[0]/2 + ((s+1)*dt+s*dt)/2 , xp, F[i,j,:]) + \
                #     w[1]*np.interp( dt*xi[0]/2 + ((s+1)*dt+s*dt)/2  , xp, F[i,j,:]) + \
                #     w[2]*np.interp( dt*xi[0]/2 + ((s+1)*dt+s*dt)/2  , xp, F[i,j,:]) \
                #     )
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

    
dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

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
tauxyk = np.zeros((mesh[0],mesh[1],nsteps+1))
tauxyp = np.zeros((mesh[0],mesh[1],nsteps+1))

step = 0 
n1,n2 = 0,0

if pbar:
    bar = progressbar.ProgressBar(max_value=nsteps)
while step <= nsteps: 
    particle_map = [[[] for _ in range(mesh[1])] for _ in range(mesh[0])]

    pos = pd.read_csv(res_dir+"position.csv."+str(step), header=None, names = ["x","y"])
    vel = pd.read_csv(res_dir+"velocity.csv."+str(step), header=None, names = ["v_x", "v_y"])
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
            for nn in range(len(particle_map[i][j])):
                # print("i = {}, j = {}, nn = {}".format(i,j,nn))
                n1 = particle_map[i][j][nn]
                ep1 = epsilon[pos_vel.loc[n1,'tipo']]
                sigma1 = sigma[pos_vel.loc[n1,'tipo']]
                rs1 = rs[pos_vel.loc[n1,'tipo']]
                m = mass[pos_vel.loc[n1,'tipo']]
                x1 = np.array([pos_vel.loc[n1,'x'], pos_vel.loc[n1,'y']])
                tauxyk[i,j,step] += -m*pos_vel.loc[n1,'v_x']*pos_vel.loc[n1,'v_y'] # falta multiplicar por 1/Volume
                # na própria célula 
                for nn2 in range(nn+1,len(particle_map[i][j])):
                    n2 = particle_map[i][j][nn2]
                    rs2 = rs[pos_vel.loc[n2,'tipo']]
                    ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                    sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                    ep = np.sqrt(ep1*ep2)
                    s = 0.5*(sigma2 + sigma1)
                    x2 = np.array([pos_vel.loc[n2,'x'], pos_vel.loc[n2,'y']])
                    r = np.sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2)
                    if r == 0:
                        print('nn {} {}'.format(nn,nn2))
                        print('x1 = {} {} | x2 = {} {}'.format(x1[0],x1[1],x2[0],x2[1]))
                        time.sleep(2)
                    coss = (x1[0]-x2[0])/r 
                    sine = (x1[1]-x2[1])/r 
                    r = r - rs1 - rs2
                    F = -(1/r**2)*(s/r)**6*(1-2*(s/r)**6)*24*ep*np.array([(x1[0]-x2[0]) - (rs1+rs2)*coss, (x1[1]-x2[1]) - (rs1+rs2)*sine])
                    tauxyp[i,j,step] += -F[1]*(x1[0] - x2[0])  # Falta multiplicar por 1/Volume
                # na celula vizinha i+1, j 
                if i+1 < region[1]:
                    print('A')
                    for nn2 in range(1,len(particle_map[i+1][j])):
                        n2 = particle_map[i+1][j][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        x2 = np.array([pos_vel.loc[n2,'x'], pos_vel.loc[n2,'y']])
                        r = np.sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2) 
                        coss = (x1[0]-x2[0])/r 
                        sine = (x1[1]-x2[1])/r 
                        r = r - rs1 - rs2
                        F = -(1/r**2)*(s/r)**6*(1-2*(s/r)**6)*24*ep*np.array([(x1[0]-x2[0]) - (rs1+rs2)*coss, (x1[1]-x2[1]) - (rs1+rs2)*sine])
                        tauxyp[i+1,j,step] += -F[1]*(x1[0] - x2[0])  # Falta multiplicar por 1/Volume
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                # na celula vizinha i+1, j+1 
                if i+1 < region[1] and j+1 < region[3]:
                    print('B')
                    for nn2 in range(1,len(particle_map[i+1][j+1])):
                        n2 = particle_map[i+1][j+1][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        x2 = np.array([pos_vel.loc[n2,'x'], pos_vel.loc[n2,'y']])
                        r = np.sqrt((x1[0] - x2[0])**2 + (x1[0] - x2[0])**2) 
                        coss = (x1[0]-x2[0])/r 
                        sine = (x1[1]-x2[1])/r 
                        r = r - rs1 - rs2
                        F = -(1/r**2)*(s/r)**6*(1-2*(s/r)**6)*24*ep*np.array([(x1[0]-x2[0]) - (rs1+rs2)*coss, (x1[1]-x2[1]) - (rs1+rs2)*sine])
                        tauxyp[i+1,j+1,step] += -F[1]*(x1[0] - x2[0])  # Falta multiplicar por 1/Volume
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                # na celula vizinha i, j+1 
                if j+1 < region[3]:
                    print('C')
                    for nn2 in range(1,len(particle_map[i][j+1])):
                        n2 = particle_map[i][j+1][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        x2 = np.array([pos_vel.loc[n2,'x'], pos_vel.loc[n2,'y']])
                        r = np.sqrt((x1[0] - x2[0])**2 + (x1[0] - x2[0])**2) 
                        coss = (x1[0]-x2[0])/r 
                        sine = (x1[1]-x2[1])/r 
                        r = r - rs1 - rs2
                        F = -(1/r**2)*(s/r)**6*(1-2*(s/r)**6)*24*ep*np.array([(x1[0]-x2[0]) - (rs1+rs2)*coss, (x1[1]-x2[1]) - (rs1+rs2)*sine])
                        tauxyp[i,j+1,step] += -F[1]*(x1[0] - x2[0])  # Falta multiplicar por 1/Volume
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                # na celula vizinha i-1, j+1 
                if i-1 >= region[0] and j+1 < region[3]:
                    print('D')
                    for nn2 in range(1,len(particle_map[i-1][j+1])):
                        n2 = particle_map[i-1][j+1][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        x2 = np.array([pos_vel.loc[n2,'x'], pos_vel.loc[n2,'y']])
                        r = np.sqrt((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2) 
                        coss = (x1[0]-x2[0])/r 
                        sine = (x1[1]-x2[1])/r 
                        r = r - rs1 - rs2
                        F = -(1/r**2)*(s/r)**6*(1-2*(s/r)**6)*24*ep*np.array([(x1[0]-x2[0]) - (rs1+rs2)*coss, (x1[1]-x2[1]) - (rs1+rs2)*sine])
                        tauxyp[i-1,j+1,step] += -F[1]*(x1[0] - x2[0])  # Falta multiplicar por 1/Volume
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
    if pbar:
        bar.update(step)
    step += 1



tauxyp = tauxyp/Vol
tauxyk = tauxyk/Vol

# vamos fazer um gráfico 

na = int(input('Enter de number of assembles (na) to average.\nThe number of steps for correlation will be nsteps/na: '))
etakk = np.zeros((mesh[0],mesh[1],int(nsteps/na)))
etakp = np.zeros((mesh[0],mesh[1],int(nsteps/na)))
etapp = np.zeros((mesh[0],mesh[1],int(nsteps/na)))
for i in range(region[0], region[1]):
    for j in range(region[2], region[3]):
        etakk[i,j,:],etakp[i,j,:],etapp[i,j,:] = greenkubo(tauxyk,tauxyp,i,j,na,dt)

ij = input('i j pro grafico: ').split()
i,j = int(ij[0]),int(ij[1])
t = np.linspace(0,t_fim/na,int(nsteps/na))
plt.figure(1)
plt.plot(t,etapp[i,j,:],label='etapp*t')
plt.plot(t,etakp[i,j,:],label='etakp*t')
plt.plot(t,etakk[i,j,:],label='etakk*t')
etatot = etakk[i,j,:]+etakp[i,j,:]+etapp[i,j,:]
plt.plot(t,etatot,'.k',label='sum')
m,b = np.polyfit(t, etatot, 1)
y = m*t + b 
plt.plot(t,y,'-k',label='Eta*t fit')
plt.legend()
plt.show()


            