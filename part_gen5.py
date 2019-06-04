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

N = 2000
L1 = 500 
L2 = 100
spac = 0.999*((L1*L2)/(N))**0.5
d1 = int(L1/spac)
d2 = int(L2/spac)
N = d1*d2
N1 = N
# Posição das partículas
p1 = np.zeros([N1,2])

#posição de referência
refe = np.array([0.01,0.01])
print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*spac+refe[0],L1, d2*spac+refe[1], L2))
print("N = {}\n".format(N))
cont = 0


for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

plt.scatter(p1[:,0],p1[:,1], label='rt1')
plt.legend()
plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))


#####################################
