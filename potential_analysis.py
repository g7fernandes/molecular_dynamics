import numpy as np 
import pandas as pd 
import os

a = input("Enter the subdomains mesh dimensions.\n")
a = a.split()
mesh = np.array([int(a[1]), int(a[2])])

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

hx = dimx/mesh[1]
hy = dimx/mesh[2]
density_map = [ [ [] ]*mesh(1) ]*mesh(2)
step = input('Choose step:\n')
while int(step) >= 0: 
    pos = pd.read_csv(res_dir+"position.csv."+step, header=None, names = ["x","y"])
    vel = pd.read_csv(res_dir+"velocity.csv."+step, header=None, names = ["v_x", "v_y"])
    n = [x for x in range(len(pos))]
    n = pd.DataFrame(n, columns=["n"]) # numero id da partícula
    pos_vel = pd.concat([n,pos,vel],axis=1)
    
    for nn in range(len(pos_vel)):
        xp = pos_vel.loc[nn,'x']//hx
        yp = pos_vel.loc[nn,'y']//hy
        if xp == mesh[1]:
            xp = xp - 1
        if yp == mesh[2]
            yp = yp - 1
        density_map[xp][yp].append(pos_vel.loc[nn,'n'])

        
    

    mean_Kenergy =  (1/2)*velocity_map**2/concentration_map
    mean_Kenergy[mean_Kenergy == np.inf] = 0
    mean_Kenergy = np.nan_to_num(mean_Kenergy)
    
    

    # gráfico 

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