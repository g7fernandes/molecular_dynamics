import numpy as np 
import os

file_name = input('Entre o arquivo:\n')
while not os.path.isfile(file_name):
    file_name = input('Entre o arquivo:\n')

data = np.genfromtxt(file_name,delimiter=',')

print('menores posições: [{}, {}] '.format( np.amin(data[:,0]), np.amin(data[:,1]) ))
print('maiores posições: [{}, {}] '.format( np.amax(data[:,0]), np.amax(data[:,1]) ))
print('medias absolutas: [{}, {}] sd: [{},{}]'.format(np.mean(abs(data[:,0])), np.mean(abs(data[:,1])), np.std(abs(data[:,0])),np.mean(abs(data[:,1])) ))
print('medias absolutas resultantes: {} sd: {}'.format( np.mean(np.sqrt(data[:,1]**2+data[:,0]**2)), np.std(np.sqrt(data[:,1]**2+data[:,0]**2)) ))

s = input('Entre os fatores de escala x e y seguidos separados por espaço:\n').split()

data[:,0] = data[:,0]*float(s[0])
data[:,1] = data[:,1]*float(s[1])

print('menores posições: [{}, {}] '.format( np.amin(data[:,0]), np.amin(data[:,1]) ))
print('maiores posições: [{}, {}] '.format( np.amax(data[:,0]), np.amax(data[:,1]) ))

le = len(data)
for i in range(le-1):
    pos1 = data[i,:]
    for j in range(i+1,le):
        pos2 = data[j,:]
        dx = abs(pos1[0]-pos2[0])
        dy = abs(pos1[1]-pos2[1])
        res = np.sqrt(dx**2+dy**2)
        if j == 1 and i == 0:
            mdx = dx
            mdy = dy 
            mres = res
        else: 
            if dx < mdx: mdx = dx
            if dy < mdy: mdy = dy
            if res < mres: mres = res

print('Menores distâncias entre dois pontos.\ndx = {}, dy = {}, absoluta = {}'.format(mdx,mdy,mres))

    
np.savetxt(file_name[0:len(file_name)-4] + '_sc.csv', data)


