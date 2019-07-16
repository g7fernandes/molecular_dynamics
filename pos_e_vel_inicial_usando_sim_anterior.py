# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
from glob import glob
import os, shutil
import configparser

dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]

config = configparser.ConfigParser()
list_files = os.listdir(res_dir)

if os.path.isfile(res_dir + 'settings.txt'):
    config.read(res_dir + 'settings.txt')
    nimpre =  int(config['global']['nimpre'].split()[0])
else:
    print('The still not converted to vtk. Probable incomplete simulation.\n'+\
          'settings.ini will be used.\n')
    config.read(dirname + '/settings.ini') 
    nimpre = len(list_files)/2    

    
N = int(config['global']['N'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
quant = []
x_files = []
v_files =[]

for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
    x_files.append(config['par_'+str(i)]['x'].split()[0])
    x_files[-1]  = x_files[-1][1:len(x_files[-1])-1] 
    v_files.append(config['par_'+str(i)]['v_file'].split()[0])
    v_files[-1]  = v_files[-1][1:len(v_files[-1])-1]
    
step = input('Extract step no. (max {}) '.format(nimpre-1))

positions = []
with open(res_dir + 'position.csv.'+step) as file:
    for line in file:
        positions.append(line)
velocities = []
with open(res_dir + 'velocity.csv.'+step) as file:
    for line in file:
        velocities.append(line)

# verificar se não é desejavel fazer backup de arquivos de posição e velocidade
a == 'n'
i = 0
for f in x_files + v_files:
    if os.path.isfile(f) and a == 'n':
        a = input('Old positions files exist. Backup? y/n')
    if a != 'n' or a != 'N':
        try:
            if i == 0:
                bkp = 'bkp'
                os.mkdir(bkp)
                i = -1
        except:
            while i >= 0:
                try: 
                    bkp = 'bkp'+str(i)
                    os.mkdir(bkp)
                    i = -1
                except:
                    i += 1
        try:                    
            shutil.move(dirname + '/' + f, dirname  + '/' + bkp + '/' + f)                
        except:
            print('file ' +f+ ' not found!')
            
k = 0            
for i in range(ntype):
    with open(x_files[i],'w') as file:
        for j in range(quant[i]):
            file.write(positions[k])
            k += 1
    
k = 0        
for i in range(ntype):
    with open(v_files[i],'w') as file:
        for j in range(quant[i]):
            file.write(velocities[k])
            k += 1        
        

