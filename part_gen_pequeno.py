# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Sun Jan 27 22:47:50 2019

# Este programa é para colocar as partículas em seus lugares iniciais

# @author: gabriel
# """

# import numpy as np
# import matplotlib.pyplot as plt
# from random import randint

# arquivo1 = 'pequeno.csv'

# N = 400
# L1 = 200
# L2 = 200
# spac = 0.999*((L1*L2)/(N))**0.5
# d1 = int(L1/spac)
# d2 = int(L2/spac)
# N = d1*d2
# N1 = N*2
# # Posição das partículas
# p1 = np.zeros([N1,2])

# #posição de referência
# refe = np.array([0.01,0.01])
# print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
# print("N = {}, spac = {}\n".format(N*2,spac))
# cont = 0


# for i in range(d1):
#     for j in range(d2):
#         p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
#         cont = cont+1
# refe = np.array([1199.9,599.9])     

# with open(arquivo1,'w') as file:
#     for i in range(N1):
#         file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

# plt.scatter(p1[:,0],p1[:,1], label='pequeno')
# plt.show()


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

Ntot = 0 
arquivo1 = 'bin_peq1.csv'

N = 1100
L1 = 400
L2 = 400
spac = 0.999*((L1*L2)/(N))**0.5
d1 = int(L1/spac)
d2 = int(L2/spac)
N = d1*d2
Ntot = N + Ntot
# Posição das partículas
p1 = np.zeros([N,2])

#posição de referência
refe = np.array([0.01,0.01])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N1 = {}, spac = {}\n".format(N,spac))
cont = 0


for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1   

# np.random.shuffle(p1)        

with open(arquivo1,'w') as file:
    for i in range(N):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

plt.scatter(p1[:,0],p1[:,1], label='pequeno1')

arquivo1 = 'bin_peq2.csv'

N = 1100
L1 = 400
L2 = 400
spac = 0.999*((L1*L2)/(N))**0.5
d1 = int(L1/spac)
d2 = int(L2/spac)
N = d1*d2
Ntot = N + Ntot
# Posição das partículas
p1 = np.zeros([N,2])

#posição de referência
refe = np.array([0.01+spac/2,0.01+spac/2])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N2 = {}, spac = {}\n".format(N,spac))
cont = 0

for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1  

with open(arquivo1,'w') as file:
    for i in range(N):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

print('Ntot = {}'.format(Ntot))
plt.scatter(p1[:,0],p1[:,1], label='pequeno2')

plt.show()