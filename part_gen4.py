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

arquivo1 = 'rt1.csv'
arquivo2 = 'rt2.csv'

N = 47704
L1 = 600 
L2 = 100
d1 = 533
d2 = 45

N1 = d1*d2
N2 = d1*d2
spac = ((600*100)/(N1+N2))**0.5
print("N1 {} N2 {}, N1 + N2 = {}".format(N1,N2,N1+N2))
# Posição das partículas
p1 = np.zeros([N1,2])
p2 = np.zeros([N2,2])



#posição de referência
refe = np.array([0.01,99.99])
cont = 0



for i in range(d1):
    for j in range(d2):
        p1[cont,:] = [i*spac+refe[0],-j*spac+refe[1]]
        cont = cont+1

refe = [0.01,00.1]   
cont = 0
for i in range(d1):
    for j in range(d2):
        p2[cont,:] = [i*spac+refe[0],j*spac+refe[1]]
        cont = cont+1

plt.scatter(p1[:,0],p1[:,1], label='rt1')
plt.scatter(p2[:,0],p2[:,1], label='rt2')
plt.legend()
plt.show()

with open(arquivo1,'w') as file:
    for i in range(N1):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

with open(arquivo2,'w') as file:
    for i in range(N2):
        file.write('{},{}\n'.format(p2[i,0],p2[i,1]))




#####################################
