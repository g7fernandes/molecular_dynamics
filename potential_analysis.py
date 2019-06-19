import numpy as np 
import pandas as pd 
import os
import configparser
import matplotlib.pyplot as plt
from glob import glob
import time

a = input("Enter the subdomains mesh dimensions.\n")
a = a.split()
mesh = np.array([int(a[0]), int(a[1])])

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
quant = []
sigma = []
epsilon = []
rs = [] # raio sólido
m = []
tipo = [0]*N

for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
    rs.append(float(config['par_'+str(i)]['rs'].split()[0]))
    sigma.append(float(config['par_'+str(i)]['sigma'].split()[0]))
    epsilon.append(float(config['par_'+str(i)]['epsilon'].split()[0]))
    m.append(float(config['par_'+str(i)]['m'].split()[0]))
j,k = 0,0
for i in range(len(quant)):
    for j in range(quant[i]):
        tipo[j+k] = i
        
    k = sum(quant[0:i+1])
tipo = pd.DataFrame(tipo, columns=["tipo"]) # numero id da partícula

hx = dimx/mesh[0]
hy = dimx/mesh[1]
particle_map = [[[] for _ in range(mesh[1])] for _ in range(mesh[0])]
density_map = np.zeros((mesh[0],mesh[1]))
kinetic_map = np.zeros((mesh[0],mesh[1]))
potential_map = np.zeros((mesh[0],mesh[1]))
step = input('Choose step:\n')
n1,n2 = 0,0
while int(step) >= 0: 
    pos = pd.read_csv(res_dir+"position.csv."+step, header=None, names = ["x","y"])
    vel = pd.read_csv(res_dir+"velocity.csv."+step, header=None, names = ["v_x", "v_y"])
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
        # print('xp = {}, yp = {}\n'.format(xp,yp))
        particle_map[xp][yp].append( pos_vel.loc[nn,'n'] )
        # print(particle_map)
        # time.sleep(0.5)
        density_map[xp,yp] += 1
        kinetic_map[xp,yp] += 0.5*np.sqrt(pos_vel.loc[nn,'v_x']**2 + pos_vel.loc[nn,'v_y']**2)*m[pos_vel.loc[nn,'tipo']]
    
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            for nn in range(len(particle_map[i][j])):
                n1 = particle_map[i][j][nn]
                ep1 = epsilon[pos_vel.loc[n1,'tipo']]
                sigma1 = sigma[pos_vel.loc[n1,'tipo']]
                rs1 = rs[pos_vel.loc[n1,'tipo']]
                # na própria célula 
                for nn2 in range(nn+1,len(particle_map[i][j])):
                    n2 = particle_map[i][j][nn2]
                    rs2 = rs[pos_vel.loc[n2,'tipo']]
                    ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                    sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                    ep = np.sqrt(ep1*ep2)
                    s = 0.5*(sigma2 + sigma1)
                    r = np.sqrt((pos_vel.loc[n1,'x'] - pos_vel.loc[n2,'x'])**2 + (pos_vel.loc[n2,'y'] - pos_vel.loc[n1,'y'])**2) - rs1 - rs2

                    potential_map[i,j] += 4*ep*((s/r)**12 - (s/r)**6)
                # na celula vizinha i+1, j 
                if i+1 < mesh[0]:
                    for nn2 in range(1,len(particle_map[i+1][j])):
                        n2 = particle_map[i+1][j][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        r = np.sqrt((pos_vel.loc[n1,'x'] - pos_vel.loc[n2,'x'])**2 + (pos_vel.loc[n2,'y'] - pos_vel.loc[n1,'y'])**2) - rs1 - rs2
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                        U = 4*ep*((s/r)**12 - (s/r)**6)
                        potential_map[i,j] += U/2
                        potential_map[i+1,j] += U/2
                # na celula vizinha i+1, j+1 
                if i+1 < mesh[0] and j+1 < mesh[1]:
                    for nn2 in range(1,len(particle_map[i+1][j+1])):
                        n2 = particle_map[i+1][j+1][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        r = np.sqrt((pos_vel.loc[n1,'x'] - pos_vel.loc[n2,'x'])**2 + (pos_vel.loc[n2,'y'] - pos_vel.loc[n1,'y'])**2) - rs1 - rs2
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                        U = 4*ep*((s/r)**12 - (s/r)**6)
                        potential_map[i,j] += U/2
                        potential_map[i+1,j+1] += U/2
                # na celula vizinha i, j+1 
                if j+1 < mesh[1]:
                    for nn2 in range(1,len(particle_map[i][j+1])):
                        n2 = particle_map[i][j+1][nn2]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        r = np.sqrt((pos_vel.loc[n1,'x'] - pos_vel.loc[n2,'x'])**2 + (pos_vel.loc[n2,'y'] - pos_vel.loc[n1,'y'])**2) - rs1 - rs2
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                        U = 4*ep*((s/r)**12 - (s/r)**6)
                        potential_map[i,j] += U/2
                        potential_map[i,j+1] += U/2
                # na celula vizinha i-1, j+1 
                if i-1 >= 0 and j+1 < mesh[1]:
                    for nn2 in range(1,len(particle_map[i-1][j+1])):
                        n2 = particle_map[i-1][j+1][nn2]
                        ep2 = epsilon[pos_vel.loc[n2,'tipo']]
                        rs2 = rs[pos_vel.loc[n2,'tipo']]
                        sigma2 = sigma[pos_vel.loc[n2,'tipo']]
                        ep = np.sqrt(ep1*ep2)
                        s = 0.5*(sigma2 + sigma1)
                        r = np.sqrt((pos_vel.loc[n1,'x'] - pos_vel.loc[n2,'x'])**2 + (pos_vel.loc[n2,'y'] - pos_vel.loc[n1,'y'])**2) - rs1 - rs2
                        if r == 0:
                            print('nn {} {}'.format(nn,nn2))
                            time.sleep(2)
                        U = 4*ep*((s/r)**12 - (s/r)**6)
                        potential_map[i,j] += U/2
                        potential_map[i-1,j+1] += U/2
                
    mean_Kenergy =  (1/2)*kinetic_map**2/density_map
    mean_Kenergy[mean_Kenergy == np.inf] = 0
    mean_Kenergy = np.nan_to_num(mean_Kenergy)
    
    X = np.linspace(0,dimx,mesh[0])
    Y = np.linspace(0,dimy,mesh[1])
    X,Y = np.meshgrid(Y,X)

    # gráfico 
    fig1 = plt.figure(1)
    kmap = plt.contourf(X,Y,mean_Kenergy,cmap="coolwarm")
    cbar1 = fig1.colorbar(kmap)
    cbar1.ax.set_ylabel('Kinetic Energy')

    fig2 = plt.figure(2)
    dmap = plt.contourf(X,Y,density_map)
    cbar2 = fig2.colorbar(dmap)
    cbar2.ax.set_ylabel('Number of density')



    # fig, ax1 = plt.subplots()
    # ax1.plot(np.linspace(0,dimx,subspaces),mean_Kenergy,'b')
    # ax1.set_xlabel("x - direction")
    # ax1.set_ylabel("Mean kinetic energy", color='b')

    # ax2 = ax1.twinx()
    # ax2.plot(np.linspace(0,dimx,subspaces),concentration_map, 'r')
    # ax2.set_ylabel("Number of density",color='r')
    # #plt.axis([0,dimx,np.min(mean_Kenergy),np.max(mean_Kenergy)])
    # fig.tight_layout()
    plt.show(block=True)

    step = input('Choose step or enter -1 to stop:\n')