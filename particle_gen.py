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

arquivo1 = 'grupo1.csv'
arquivo2 = 'grupo2.csv'

N1 = 1600
N2 = 6400
# Posição das partículas
p1 = np.zeros([N1,2])
p2 = np.zeros([N2,2])


d1 = 160
d2 = 40


spac = 2**(1/6)

#posição de referência
refe = np.array([97,100])
cont = 0



for i in range(d2):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

refe = [30,50]   
cont = 0
for i in range(d1):
    for j in range(d2):
        p2[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

plt.scatter(p1[:,0],p1[:,1], label='grupo1')
plt.scatter(p2[:,0],p2[:,1], label='grupo2')
plt.legend()
plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

with open(arquivo2,'w') as file:
    for i in range(N2):
        file.write('{},{}\n'.format(p2[i,0],p2[i,1]))




#####################################
