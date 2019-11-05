import numpy as np 
import os
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)


lab = 'arq'
aux = 0
outro = 'a'
while outro != 'n':
    file_name = input('Entre o arquivo:\n')
    while not os.path.isfile(file_name):
        file_name = input('Entre o arquivo:\n')

    data = np.genfromtxt(file_name,delimiter=',')
    ax.scatter(data[:,0],data[:,1], label=lab+str(aux))
    
    outro = input('adicionar outro arquivo? s/n ')
    aux += 1
    
ax.set_aspect(aspect=1)

plt.legend()
plt.show()


    
