import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

sigma_ref = 0.341

delta = 0.005
sig0 = 1
ep0 = 1
rs = 5 # raio solido
rm = .245 / sigma_ref #raio do atomo
# .245 é a distância entra atomos em estrutura hexagonal rombica

x = y = np.arange(-(3*sig0+rs), 3*sig0+rs, delta)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)
np.place(r,r<rs,np.nan)
gamma = np.arctan(X/Y)


A = (0.071/sigma_ref)/2 # no epsilon
B = (0.071/sigma_ref)/2 # no sigma
alpha = beta = 2*np.pi*rs/rm
ph = 2*np.pi/alpha
#beta = 6

V = 4*ep0*(1+A*np.sin(alpha*gamma+ph))* ( (sig0*(1+B*np.sin(beta*gamma))/(r-rs) )**12 - (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6  )

Fr = 24*ep0*(1+A*np.sin(alpha*gamma+ph))*(1/(r-rs)**2) *  ( (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 * (1-2* (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 )) *  (r-rs)

Fg = (1/(r-rs))* ( (4*ep0*(A*np.sin(alpha*gamma+ph) + 1 ))* (sig0/(r-rs))**6*(B*np.sin(beta*gamma)+1)**5 * \
     ( (sig0/(r-rs))**6 * 12*beta*B*np.cos(beta*gamma)*(B*np.sin(beta*gamma)+1)**6 - 6*beta*B*np.cos(beta*gamma)) + \
    4*alpha*A*ep0*np.cos(alpha*gamma+ph)*(sig0/(r-rs))**6*(B*np.sin(beta*gamma)+1)**6 * \
    ((sig0/(r-rs))**6 * (B*np.sin(beta*gamma)+1)**6 - 1) ) 


# Fg = (1/r)* ( (4*ep0*(A*np.sin(alpha*gamma) + 1 ))* \
#      ( (sig0/r)**12 * 12*beta*B*np.cos(beta*gamma)*(B*np.sin(beta*gamma)+1)**11 - (sig0/r)**6 * 6*beta*B*np.cos(beta*gamma)*(B*np.sin(beta*gamma)+1)**5) + \
#     4*alpha*A*ep0*np.cos(alpha*gamma)*((sig0/r)**12 * (B*np.sin(beta*gamma)+1) - (sig0/r)**6*(B*np.sin(beta*gamma)+1)**6) ) 



# pot 
np.place(V,abs(V) > 20, np.nan)

plt.figure()
plt.plot(X[int(len(X)/2),:] ,V[int(len(X)/2),:])
plt.title("potential")


fig1, ax1 = plt.subplots()
cs1 = ax1.contourf(X, Y, V,levels=20)
cbar = fig1.colorbar(cs1)
ax1.set_title("potential")

#force 

np.place(Fr,abs(Fr) > 5000, np.nan)

# plt.figure()
# plt.plot(X[int(len(X)/2),:] ,Fr[int(len(X)/2),:])
# plt.title("force")


cmap = plt.get_cmap('PiYG')
fig2, ax2 = plt.subplots()
cs2 = ax2.contourf(X, Y, Fr,levels=20, cmap=cmap)
cbar = fig2.colorbar(cs2)
ax2.set_title('radial force')


np.place(Fg,abs(Fg) > 5000, np.nan)
fig3, ax3 = plt.subplots()
cs3 = ax3.contourf(X, Y, Fg, levels=20, cmap=cmap)
cbar = fig3.colorbar(cs3)
ax3.set_title('tangential force')


plt.show()