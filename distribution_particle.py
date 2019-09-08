import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import configparser

def trapz(x,y):
    res = 0
    for i in range(1,len(x)):
        res += (y[i-1]+y[i])*(x[i-1]+x[i])/2
    return res

# LÊ A CONFIGURAÇÃO

config = configparser.ConfigParser()
config.read('settings.ini')
print('Reading settings.ini')
N = int(config['global']['N'].split()[0])
nimpre =  int(config['global']['nimpre'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
nimpre_init = int(config['global']['nimpre_init'].split()[0])
dimx = int(config['global']['dimX'].split()[0])
dimy = int(config['global']['dimY'].split()[0])

quant = []


for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))

print("Initial step: {}. Final step: {}. Number of particle groups (ntype): {}\n".format(nimpre_init,nimpre,ntype))

pgroup = input("Enter the particle group (>= 0, < ntype): ")
stepini = input("Enter the initial step: ")
stepfim = input("Enter the final step: ")

x = np.zeros((stepfim - stepini, quant[int(pgroup)], 2)

# LÊ OS VTK

grupo = "./result/grupo" + pgroup + "_"

reader = vtk.vtkXMLUnstructuredGridReader()
nfiles = 0
for file in range(stepini,stepfim)
    reader.SetFileName(grupo+str(file))
    reader.Update()
    data = reader.GetOutput()

    points = data.GetPoints()
    x[nfiles,:,:] = vtk_to_numpy(points.GetData())[:,0:2]# coluna 1 = x, coluna 2  = y
    nfiles += 1
    #u = vtk_to_numpy(data.GetPointData().GetArray(0)) # colunas: Vx, Vy, Vz

# Todas as posições carregadas

aux1 = input("Enter the number of subdivisions and the direction (x or y).\nIf direction not given, will be x.\n")).split()
direc = 0
if len(aux1) > 1:
    if aux[0] == 'y':
        direc = 1


ndensiy = np.zeros((nfiles,aux1))
dhist = dimx/int(aux1[0])  # distancia no histograma

for file in range(nfiles):
    for i in range(quant[int(pgroup)]):
        ndensiy[file, (x[file,i,direc]//dhist) ] += 1   


n2density = (np.sum(ndensity,axis=0)/nfiles)/(quant[int(pgroup)])
xx = np.linspace(0,dimx, int(aux1[0]))
plt.plot(n2density)

du = quant[int(pgroup)]/dimx #distribuição uniforme


DS = trapz(x,n2density*np.log(n2density/du))





#N = reader.GetNumberOfPoints() # numero de pontos
#npa = reader.GetNumberOfPointArrays() # numero de vetores de ponto

# n_arrays = reader.GetNumberOfPointArrays()
# print("Available data:\n")
# for i in range(n_arrays):
#     print(reader.GetPointArrayName(i))



