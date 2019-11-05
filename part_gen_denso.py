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

#### sort the particles into regions ####
def density_map(x,dimx,dimy,tam):

    xp = x[:,0] // tam # This gives the location of each particle
    yp = x[:,1] // tam
    dmap = [[[] for _ in range(int(dimy/tam)+1)] for _ in range(int(dimx/tam+1))]
    for i in range(len(x)):
        dmap[int(xp[i])][int(yp[i])].append(i)
    return dmap
#########################################

tol = 1

ff = 1 # formato da região ff = x/y
Lx = 100
arquivo1 = 'molp.csv'
arquivo2 = 'parp.csv'

N1 = 1000
N2 = 10

L1 = Lx
L2 = Lx/ff
print("Formato da região: {} x {}".format(L1,L2))
spac = 2**(1/6)
N1 = int((int((tol*L1-spac)/spac)+ (int((tol*L1-spac)/spac)-2))*int((tol*L2-spac)/(spac**2*(3/4)))/2)
aux = '1'
while aux != '':
    print("Espaçamento com N1 = {} partículas: {}.\n".format(N1,spac))
    aux = input("Entre outro valor de distância entre particulas se desejar: ".format(2**(1*6)))
    if aux != '':
        spac = eval(aux)
        N1 = int((int(tol*L1/spac)+ (int(tol*L1/spac)-2))*int(tol*L2/(spac**2*(3/4)))/2)
    print("Espaçamento com N1 = {} partículas: {}.\n".format(N1,spac))

# spac = 1
p1 = []

#posição de referência
refe = np.array([spac/2,spac/2])

cont = 0
x,y = refe[0], refe[1]
flag = True
while y < L2*tol:
    while x < L1*tol:
        p1.append([x,y])
        x += spac
        cont += 1
    y += spac*3/4
    if flag:
        x = refe[0] + spac/2
        flag = False
    else:
        x = refe[0]
        flag = True

p1 = np.array(p1)


r2 = float(input(">> Raio da partícula 2: "))

if N2 > 0:
    d1 = int(np.sqrt(N2*(L1*tol-2*r2)/(L2*tol-2*r2)))
    d2 = int(N2/d1)
    sx = (L1*tol-2*r2)/d1
    sy = (L2*tol-2*r2)/d2

    N2 = d1*d2
    # Posição das partículas
    p2 = np.zeros([N2,2])

    #posição de referência
    refe = np.array([0.0005+r2, 0.0005+r2])
    print("d1*spac = {}, L1 = {}, d2*spac = {}, L2 = {}\n".format(d1*sx+refe[0],L1, d2*sy-refe[1], L2))
    print("N2 = {}, spac = {} {}\n".format(N2,sx,sy))
    cont = 0

    for i in range(d1):
        for j in range(d2,0,-1):
            p2[cont,:] = [i*sx+refe[0],j*sy+refe[1]]
            cont = cont+1 

    cx = (np.min(p2[:,0]) + np.max(p2[:,0]))/2
    cy = (np.min(p2[:,1]) + np.max(p2[:,1]))/2
    p2[:,0] = p2[:,0] + L1/2 - cx 
    p2[:,1] = p2[:,1] + L2/2 - cy

    dmap = density_map(p1,L1,L2,r2)

    for k in range(len(p2)):
        m = [int(p2[k,0]//r2), int(p2[k,1]//r2)]
        for i in range(m[0]-1,m[0]+2):
            for j in range(m[1]-1,m[1]+2):
                if (i > -1 and j > -1 and i <= int(L1/r2) and j <= int(L2/r2)):
                    for n in range(len(dmap[i][j])):
                        x1 = p1[dmap[i][j][n],:]
                        r = np.sqrt((x1[0]-p2[k,0])**2 + (x1[1]-p2[k,1])**2)
                        if r <= r2:
                            p1[dmap[i][j][n],:] = np.array([np.nan, np.nan])

    p1 = p1[~np.isnan(p1[:,0]),:]

# Plota e salva arquivos
fig = plt.figure()
ax = fig.add_subplot(111)

np.random.shuffle(p1)        

with open(arquivo1,'w') as file:
    for i in range(len(p1)):
        file.write('{},{}\n'.format(p1[i,0],p1[i,1]))

ax.scatter(p1[:,0],p1[:,1], s=(np.pi*1**2),label='arquivo1')

if N2 > 0:
    np.random.shuffle(p2)        

    with open(arquivo2,'w') as file:
        for i in range(len(p2)):
            file.write('{},{}\n'.format(p2[i,0],p2[i,1]))

    Ntot = N1 + N2
    ax.scatter(p2[:,0],p2[:,1], alpha = 0.5, s=(np.pi*r2**2),label='arquivo2')

    print('N1 = {}, N2 = {}, Ntot = {}\n'.format(len(p1),len(p2),len(p1)+len(p2)))
    print("Arquivos {} e {}\n".format(arquivo1,arquivo2))
else:
    print('N1 = {}\n'.format(len(p1)))
    print("Arquivos {} e {}\n".format(arquivo1,arquivo2)) 

ax.set_aspect(aspect=1)
ax.legend()
plt.show()



