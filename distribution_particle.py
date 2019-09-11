import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import configparser
import pandas as pd
import os.path
from glob import glob 
import matplotlib.pyplot as plt
try:
    import vtk
except:
    print("Vtk module not found! If using anaconda try:\nconda install -c anaconda vtk\n OR\n pip install vtk")
## --------  Funções --------- ## 
def trapz(x,y):
    res = 0
    for i in range(1,len(x)):
        res += (y[i-1]+y[i])*(x[i-1]+x[i])/2
    return res

#---------------------------------------------------------------#

# escolhe a pasta

dirname = os.getcwd() #os.path.dirname(os.path.abspath(__file__))
dirlist = glob(dirname + "/*/")
print("Choose a folder there the results are contained:\nNo | Folder")
for a in range(len(dirlist)):
    print("{} | {}\n".format(a,dirlist[a]))
a = int(input("Enter the number of the folder\n"))
res_dir = dirlist[a]
res_name = res_dir.split('/')[-2] # nome da pasta com o resultado

# LÊ A CONFIGURAÇÃO

config = configparser.ConfigParser()
config.read(res_dir + 'settings.txt')
print('Reading settings.txt')
N = int(config['global']['N'].split()[0])
nimpre =  int(config['global']['nimpre'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
nimpre_init = int(config['global']['nimpre_init'].split()[0])
dimx = float(config['global']['dimX'].split()[0])
dimy = float(config['global']['dimY'].split()[0])

quant = []


for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))

print("Initial step: {}. Final step: {}. Number of particle groups (ntype): {}\n".format(nimpre_init,nimpre,ntype))

pgroup = input("Enter the particle group (>= 0, < ntype):[{}] ".format(ntype-1))
if (pgroup == ''):
    pgroup = str(ntype-1)
    print("pgroup = {}\n".format(pgroup))
stepini = input("Enter the initial step:[0] ")
if (stepini == ''):
    stepini = 0
    print("stepini = 0\n") 
else:
    stepini = int(stepini)
stepfim = input("Enter the final step :[{}] ".format(nimpre))
if (stepfim == ''):
    stepfim = nimpre
    print("stepfim = {}".format(nimpre))
else:
    stepfim = int(stepfim)

aux1 = input("Number of files to use (will skip some if less than given above):[{}] ".format(nimpre/2))
if aux1 == '':
    aux1 = nimpre/2
    print("jump = {}\n".format(int(nimpre/2)))
else:
    aux1 = int(aux1)    
jump = int((stepfim - stepini)/aux1)

x = np.zeros((stepfim - stepini, quant[int(pgroup)], 2))

# LÊ OS VTK

grupo = res_dir + "grupo" + pgroup + "_"

reader = vtk.vtkXMLUnstructuredGridReader()
nfiles = 0
proximo = stepini 
for file in range(stepini,stepfim):
    if file == proximo:
        reader.SetFileName(grupo+str(file)+'.vtu')
        reader.Update()
        data = reader.GetOutput()

        points = data.GetPoints()
        x[nfiles,:,:] = vtk_to_numpy(points.GetData())[:,0:2]# coluna 1 = x, coluna 2  = y
        nfiles += 1
        proximo += jump
        
    #u = vtk_to_numpy(data.GetPointData().GetArray(0)) # colunas: Vx, Vy, Vz

# Todas as posições carregadas

aux1 = input("Enter the number of subdivisions and the direction (x or y).\nIf direction not given, will be x.\n [10 x]").split()
if len(aux1) == 0:
    aux1 = ['10']
    print("There will be 10 subdivisions.\n")
direc = 0
if len(aux1) > 1:
    if aux1[0] == 'y':
        direc = 1


ndensiy = np.zeros((nfiles,int(aux1[0])))
dhist = dimx/int(aux1[0])  # distancia no histograma

for file in range(nfiles):
    for i in range(quant[int(pgroup)]):
        ndensiy[file, int(x[file,i,direc]//dhist) ] += 1   


ndensity2 = (np.sum(ndensiy,axis=0)/nfiles)/(quant[int(pgroup)])
xx = np.linspace(0,dimx, int(aux1[0]))
plt.plot(ndensity2)

du = quant[int(pgroup)]/dimx #distribuição uniforme


DS = trapz(xx,ndensity2*np.log(ndensity2/du))
print("Relative distribution relative to the Uniform: {}\n".format(DS))
plt.show()


#N = reader.GetNumberOfPoints() # numero de pontos
#npa = reader.GetNumberOfPointArrays() # numero de vetores de ponto

# n_arrays = reader.GetNumberOfPointArrays()
# print("Available data:\n")
# for i in range(n_arrays):
#     print(reader.GetPointArrayName(i))



##-- SALVA dados em arquivos csv --##
if os.path.exists('Distributions.csv'):
    csv_input = pd.read_csv('Distributions.csv')
    headers_list = list(csv_input)

    laux1 = True # variavel auxiliar lógica
    while laux1:
        for name in headers_list:
            if name == res_name:
                print("A file with name {} already exists.")
                res_name = input('Enter a new name or ctrl+C to abort:\n')
                laux1 = True 
                break
            else:
                laux1 = False
    
    csv_input[res_name] = ndensity2
else:
    csv_input = pd.DataFrame(ndensity2,columns=[res_name])


csv_input.to_csv('Distributions.csv', index=False)

#--#


f=open('Distrib_entropies.csv', "a+")
f.write("{}, {}\n".format(res_name,DS))





