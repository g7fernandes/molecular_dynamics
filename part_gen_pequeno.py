# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Sun Jan 27 22:47:50 2019

# Este programa é para colocar as partículas em seus lugares iniciais

# @author: gabriel
# """



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 22:47:50 2019

Este programa é para colocar as partículas em seus lugares iniciais

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
from random import randint

ff = 1 # formato da região ff = x/y
Lx = 400
arquivo2 = 'poucas1.csv'
arquivo1 = 'poucas2.csv'

ratio = 0
N1 = 500
N2 = 100

L1 = Lx
L2 = Lx/ff
spac = 0.999*((L1*L2)/(N1))**0.5


d1 = int(L1/spac)
d2 = int(L2/spac)
N1 = d1*d2

sx = L1*.99/(d1+1) 
sy = L2*.99*ratio/(d2+1)
# Posição das partículas
p1 = np.zeros([N1,2])

#posição de referência
refe = np.array([0.01,0.01])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*sx+refe[0],L1, d2*sy+refe[1], L2))
print("N1 = {}, spac = {}, {}\n".format(N1,sx,sy))
cont = 0


for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*sx+refe[0],j*sy+refe[1]]
        cont = cont+1   

np.random.shuffle(p1)        

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

plt.scatter(p1[:,0],p1[:,1], label='pequeno1')


L1 = Lx
L2 = Lx/ff
spac = 0.99*((L1*L2)/(N2))**0.5


d1 = int(L1/spac)
d2 = int(L2/(spac))

sx = L1*.99/(d1+1) 
sy = L2*.99*(1-ratio)/(d2+1)

N2 = d1*d2
# Posição das partículas
p1 = np.zeros([N2,2])

#posição de referência
refe = np.array([0.01+sx/2, L2 + 0.01+sy/2])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*sx+refe[0],L1, d2*sy-refe[1], L2))
print("N2 = {}, spac = {} {}\n".format(N2,sx,sy))
cont = 0

for i in range(d1):
    for j in range(d2,0,-1):
        p1[cont,:] = [i*sx+refe[0],-j*sy+refe[1]]
        cont = cont+1  

np.random.shuffle(p1)        

with open(arquivo2,'w') as file:
    for i in range(N2):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

Ntot = N1 + N2
print('Ntot = {}'.format(Ntot))
plt.scatter(p1[:,0],p1[:,1], label='pequeno2')

plt.legend()
plt.show()



