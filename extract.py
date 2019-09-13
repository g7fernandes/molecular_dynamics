# Extract the particles that are in a region

import numpy as np 
import matplotlib.pyplot as plt 
from vtk.util.numpy_support import vtk_to_numpy
try:
    import vtk
except:
    print("Vtk module not found! If using anaconda try:\nconda install -c anaconda vtk\n OR\n pip install vtk")

# Lê texto csv
#xfile = input("Enter position file:\n")
#vfile = input("Enter velocity file:\n")
# Lê vtl
vtkfile = input("Enter .vtu file:\n")
limits = input("Enter limits to extract x_min x_max y_min y_max:\n").split()
limits = [float(l) for l in limits]

# Lê arquivo do VTK
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(vtkfile)
reader.Update()
data = reader.GetOutput()
points = data.GetPoints()
x = vtk_to_numpy(points.GetData())[:,0:2]
v = vtk_to_numpy(data.GetPointData().GetArray(0))[0,0:2] # colunas: Vx, Vy, Vz

# Estes são para casos em que desejar ler arquivo txt
#x = np.genfromtxt(xfile, delimiter=',')
#v = np.genfromtxt(vfile, delimiter=',')

x_out = [] 
v_out = []
for i in len(x):
    if x[i,0] > limits[0] and x[i,0] < limits[1] and x[i,1] > limits[2] and x[i,1] < limits[3]:
        x_out.append([x[i,0], x[i,1]])
        v_out.append([v[i,0], v[i,1]])

x_out = np.array(x_out)
v_out = np.array(v_out)

print("Number of particles: {}\n".format(len(x_out)))
print("Mean velocity: {}\n".format(np.mean(x_out[:,0]**2 + x_out[:,1]**2))
ax.scatter(x_out[:,0],x_out[:,1])
ax.set_aspect(aspect=1)
ax.set_xlim(limits[0],limits[1])ax.set_xlim(limits[0],limits[1])
ax.set_ylim(limits[2],limits[3])


ax.set_xlim(limits[0],limits[1])
ax.set_ylim(limits[2],limits[3])

plt.show()
sug = vtkfile.split('/')[0]

fname = input("Enter file names to save without .csv. [{}]\n".format(sug))
if fname == '':
    fname = sug 

np.savetxt('x_' + fname + '.csv', x_out)
np.savetxt('v_' + fname + '.csv', v_out)





