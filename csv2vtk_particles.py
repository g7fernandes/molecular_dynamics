#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 13:03:06 2019

Este programa lê arquivos CSV onde cada arquivo é um passo de tempo e cada 
linha é uma partícula e exporta em vtk o reultado. 

@author: gabriel
"""

from evtk.hl import pointsToVTK
import numpy as np
import csv, os, shutil
import configparser

via = os.getcwd()

aux = True
folder = 'result'
i = 1

while aux:
    try:
        os.mkdir(folder)
        aux = False
    except FileExistsError:
        print('Folder {} exists'.format(folder))
        opt = input('Overwrite? [y/n]')
        if opt == 'y' or 'opt' == 'Y':
            aux =False 
        else:
            folder = input('Enter new folder name:\n')
            if len(folder) == 0:  
                folder = 'result' + '_' + str(i) 
                aux = True
                i = i+1       

print('The results will be saved at {}'.format(folder))

#Ler a quantidade de arquivos de settings.

config = configparser.ConfigParser()
config.read('settings.ini')

N = int(config['global']['N'].split()[0])
nimpre =  int(config['global']['nimpre'].split()[0])
ntype = int(config['global']['Ntype'].split()[0])
quant = []
rs = [] # raio sólido

for i in range(ntype):
    quant.append(int(config['par_'+str(i)]['quantidade'].split()[0]))
    rs.append(float(config['par_'+str(i)]['rs'].split()[0]))

a = os.listdir('temp')
if len(a)/2 < nimpre:
    print("Propable incomplete execution. Processing {} files\n.".format(len(a)/2))
    nimpre = int(len(a)/2)-1

tipo = np.zeros(N)
rsol = np.zeros(N)
j,k = 0,0
for i in range(len(quant)):
    for j in range(quant[i]):
        tipo[j+k] = i
        rsol[j+k] = rs[i]
    k = quant[i]
        
# for i in range(len(rsol)):
#     print("{} {}".format(i,rsol[i]))

# Ler os arquivos e colocá-los em vetores

x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)

vx = np.zeros(N)
vy = np.zeros(N)
vz = np.zeros(N)

cx = np.zeros(N)
cy = np.zeros(N)
cz = np.zeros(N)

try:
    shutil.move(via+'/settings.txt',via+'/'+folder+'/settings.txt') 
except:
    print('No settings file found!\n')

for fnum in range(0,nimpre+1):
    with open('temp/position.csv.'+str(fnum),encoding='utf-8') as file_locus:
        csv_lector = csv.reader(file_locus,delimiter = ',')
        i = 0
        for linea in csv_lector:
            x[i] = linea[0]
            y[i] = linea[1]
            i = i+1
    shutil.move(via+'/temp/position.csv.'+str(fnum),via+'/'+folder+'/position.csv.'+str(fnum)) 
    
    with open('temp/velocity.csv.'+str(fnum),encoding='utf-8') as file_velocitas:
        csv_lector = csv.reader(file_velocitas,delimiter = ',')
        i = 0
        for linea in csv_lector:
            vx[i] = linea[0]
            vy[i] = linea[1]
            
            i = i+1
            
    shutil.move(via+'/temp/velocity.csv.'+str(fnum),via+'/'+folder+'/velocity.csv.'+str(fnum))
    
    #with open('cell.csv.'+str(fnum),encoding='utf-8') as file_velocitas:
        #csv_lector = csv.reader(file_velocitas,delimiter = ',')
        #i = 0
        #for linea in csv_lector:
            #cx[i] = int(float(linea[0]))
            #cy[i] = int(float(linea[1]))
            
            #i = i+1
            
    #shutil.move(via+'/cell.csv.'+str(fnum),via+'/'+folder+'/cell.csv.'+str(fnum))    
    
    pointsToVTK(via+'/'+folder +'/result_'+str(fnum), x, y, z, data = {"Vx" : vx, "Vy" : vy, "Tipo" : tipo, "raio_solido" : rsol })        #, "Cellx" : cx, "Celly" : cy

#shutil.rmtree('temp')
        
        
