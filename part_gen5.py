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

arquivo1 = 'grupo0.csv'

N = 1000
L1 = 200
L2 = 600
spac = 0.999*((L1*L2)/(N))**0.5
d1 = int(L1/spac)
d2 = int(L2/spac)
N = d1*d2
N1 = N*2
# Posição das partículas
p1 = np.zeros([N1,2])

#posição de referência
refe = np.array([0.01,0.01])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N = {}, spac = {}\n".format(N*2,spac))
cont = 0


for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1
refe = np.array([1199.9,599.9])        
for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [-i*spac+refe[0],-j*spac+refe[1]]
        cont = cont+1

np.random.shuffle(p1)

plt.scatter(p1[:,0],p1[:,1], label='grupo0')
# plt.legend()
# plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))


#####################################

arquivo1 = 'p_g.csv'

N = 20
L1 = 1200
L2 = 520
spac = 1.1*((L2)/(N))
d1 = 1 #int(L1/spac)
d2 = N #int(L2/spac)
N = d1*d2
N1 = N
# Posição das partículas
p1 = np.zeros([N1,2])

#posição de referência
refe = np.array([600,20])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N = {}, spac= {}\n".format(N,spac))
cont = 0


for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

plt.scatter(p1[:,0],p1[:,1], label='pg')
#plt.legend()
#plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

#######################

arquivo1 = 'p_p.csv'

N = 20
L1 = 1200
L2 = 520
spac = 1.1*((L2)/(N))
d1 = 1 #int(L1/spac)
d2 = N #int(L2/spac)
N = d1*d2
N1 = N
# Posição das partículas
p2 = np.zeros([N1,2])

#posição de referência
refe = np.array([600,20-spac/2])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N = {}, spac= {}\n".format(N,spac))
cont = 0


for i in range(d1):
    for j in range(d2):
        p2[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

plt.scatter(p2[:,0],p2[:,1], label='pp')
plt.legend()
plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p2[i,0],p2[i,1]))
