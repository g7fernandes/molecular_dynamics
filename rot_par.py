import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
 

def abssin(x):
    return np.sin(x)/abs(np.sin(x))

maxF = 40
sigma_ref = 0.341

delta = 0.005
sig0 = 1
ep0 = 4
rs = 5 # raio solido
rm = .245 / sigma_ref #raio do atomo
# .245 é a distância entra atomos em estrutura hexagonal rombica

x = y = np.arange(0, 3*sig0+rs, delta)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)
np.place(r, r < rs, np.nan)
gamma = np.arctan2(Y,X)


A = 0.0001 #ep0/2 #(0.071/sigma_ref)/4  # no epsilon
B = 0.0001 #(0.071/sigma_ref)/8 # no sigma
alpha = beta = (2*np.pi*rs/rm)
ph = (2*np.pi/alpha)
print("A {}, B {}, alpha {}, fase {}, periodo {}\n".format(A, B, alpha ,ph,2*np.pi/alpha))

print("pr = {} {} {} {} {}\n".format(A,B,alpha,beta,ph))
rmolec = rs 
gmolec = np.linspace(0,2*np.pi/4,round((2*np.pi*rs/rm)/4))
#beta = 6

V = 4*ep0*(1+A*np.sin(alpha*gamma+ph))* ( (sig0*(1+B*np.sin(beta*gamma))/(r-rs) )**12 - (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6  )
#V = 4*ep0*(1+A*np.sin(alpha*gamma+ph))* ( (sig0*(1+B*abs(np.sin(beta*gamma)))/(r-rs) )**12 - (sig0*(1+B*abs(np.sin(beta*gamma)))/(r-rs))**6  )

#V = 4*ep0*(1+A*np.sin(alpha*gamma+ph))* ( (sig0/(r-rs*(1+B*np.sin(beta*gamma))) )**12 - (sig0/(r-rs*(1+B*np.sin(beta*gamma))))**6  )

#Fr = 24*ep0*(1+A*np.sin(alpha*gamma+ph))*(1/(r-rs)**2) *  ( (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 * (1-2* (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 )) *  (r-rs)
# Fg = (1/(r-rs))* ( (4*ep0*(A*np.sin(alpha*gamma+ph) + 1 ))* (sig0/(r-rs))**6*(B*(np.sin(beta*gamma))+1)**5 * \
#      ( (sig0/(r-rs))**6 * 12*beta*B*np.cos(beta*gamma)*(B*(np.sin(beta*gamma))+1)**6 - 6*beta*B*np.cos(beta*gamma)) + \
#     4*alpha*A*ep0*np.cos(alpha*gamma+ph)*(sig0/(r-rs))**6*(B*(np.sin(beta*gamma))+1)**6 * \
#     ((sig0/(r-rs))**6 * (B*(np.sin(beta*gamma))+1)**6 - 1) ) 

Fr = 24*ep0*(1+A*sin(alpha*gamma+ph))*(1/(r-rs)**2) * ((sig0*(1+B*sin(beta*gamma))/(r-rs))**6 *  \
                                    (1-2* (sig0*(1+B*sin(beta*gamma))/(r-rs))**6 )) * (r-rs)

Fg = (1/(r-rs))* ( (4*ep0*(A*np.sin(alpha*gamma+ph) + 1 ))* (sig0/(r-rs))**6*(B*np.sin(beta*gamma)+1)**5 * \
     ( (sig0/(r-rs))**6 * 12*beta*B*np.cos(beta*gamma)*(B*np.sin(beta*gamma)+1)**6 - 6*beta*B*np.cos(beta*gamma)) + \
    4*alpha*A*ep0*np.cos(alpha*gamma+ph)*(sig0/(r-rs))**6*(B*np.sin(beta*gamma)+1)**6 * \
    ((sig0/(r-rs))**6 * (B*np.sin(beta*gamma)+1)**6 - 1) ) 

# Fr = 24*ep0*(1+A*np.sin(alpha*gamma+ph))*(1/(r-rs)**2) *  ( (sig0*(1+B*abs(np.sin(beta*gamma)))/(r-rs))**6 * (1-2* (sig0*(1+B*abs(np.sin(beta*gamma)))/(r-rs))**6 )) *  (r-rs)

# Fg = (1/(r-rs))* ( (4*ep0*(A*np.sin(alpha*gamma+ph) + 1 ))* (sig0/(r-rs))**6*(B*abs(np.sin(beta*gamma))+1)**5 * \
#      ( (sig0/(r-rs))**6 * 12*beta*B*np.cos(beta*gamma)*abssin(beta*gamma)*(B*abs(np.sin(beta*gamma))+1)**6 - 6*beta*B*np.cos(beta*gamma)*abssin(beta*gamma)) + \
#     4*alpha*A*ep0*np.cos(alpha*gamma+ph)*(sig0/(r-rs))**6*(B*abs(np.sin(beta*gamma))+1)**6 * \
#     ((sig0/(r-rs))**6 * (B*abs(np.sin(beta*gamma))+1)**6 - 1) ) 

# pot 
np.place(V,abs(V) > 20, np.nan)
np.place(Fr,abs(Fr) > maxF, np.nan)
np.place(Fg,abs(Fg) > maxF, np.nan)

plt.figure()
plt.plot(X[int(len(X)/2),:] ,V[int(len(X)/2),:])
plt.title("potential")


fig1, ax1 = plt.subplots()
cs1 = ax1.contourf(X, Y, V,levels=20)
cbar = fig1.colorbar(cs1)
ax1.set_title("potential")
ax1.plot(rmolec*np.sin(gmolec), rmolec*np.cos(gmolec),'ko')

#force 



# plt.figure()
# plt.plot(X[int(len(X)/2),:] ,Fr[int(len(X)/2),:])
# plt.title("force")


cmap = plt.get_cmap('PiYG')
fig2, ax2 = plt.subplots()
cs2 = ax2.contour(X, Y, Fr,levels=20, colors='k' )

# ax2.set_title('radial force')



#fig3, ax3 = plt.subplots()
cs3 = ax2.contourf(X, Y, Fg, levels=20, cmap=cmap)
#cbar = fig3.colorbar(cs3)

cbar1 = fig2.colorbar(cs3)

ax2.plot(rmolec*np.sin(gmolec), rmolec*np.cos(gmolec),'ko')
ax2.set_title('tangential force')

#---------------------------------------------------------#

fig3, ax3 = plt.subplots()
cs3 = ax3.contour(X, Y, Fg, levels=20, colors='k')

cmap = plt.get_cmap('PiYG')

cs2 = ax3.contourf(X, Y, Fr,levels=20, cmap=cmap)
cbar2 = fig3.colorbar(cs2)
ax3.plot(rmolec*np.sin(gmolec), rmolec*np.cos(gmolec),'ko')
ax3.set_title('radial force')




# --------------------------------------------------------#
delta = 0.05
x = np.arange(-1.5, 1.5, delta)
y = np.arange(5.50, 8, delta)
# x = y = np.arange(0, 3*sig0+rs, delta)
X, Y = np.meshgrid(x, y)

r = np.sqrt(X**2 + Y**2)
np.place(r, r < rs, np.nan)
gamma = np.arctan2(Y,X)

Fr = 24*ep0*(1+A*np.sin(alpha*gamma+ph))*(1/(r-rs)**2) *  ( (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 * (1-2* (sig0*(1+B*np.sin(beta*gamma))/(r-rs))**6 )) *  (r-rs)

Fg = (1/r)* ( (4*ep0*(A*sin(alpha*gamma+ph) + 1 ))* (sig0/r)**6*(B*sin(beta*gamma)+1)**5 * \
                                    ((sig0/r)**6 * 12*beta*B*cos(beta*gamma)*(B*sin(beta*gamma)+1)**6 - 6*beta*B*cos(beta*gamma)) + \
                                    4*alpha*A*ep0*cos(alpha*gamma+ph)*(sig0/r)**6*(B*sin(beta*gamma)+1)**6 * \
                                    ((sig0/r)**6 * (B*sin(beta*gamma)+1)**6 - 1)) 
np.place(Fr,abs(Fr) > maxF, np.nan)
np.place(Fg,abs(Fg) > maxF, np.nan)

fig4, ax4 = plt.subplots()
ax4.set_title("Force")
Fy =  (-Fr)*np.sin(gamma) + -Fg*np.cos(-gamma)
Fx =  (-Fr)*np.cos(gamma) + -Fg*np.sin(gamma)
M = np.hypot(Fx, Fy)
Q = ax4.quiver(X, Y, Fx, Fy, M) #  units='x', pivot='tip', width=0.022, scale=1 / 0.15)
cb = fig4.colorbar(Q)
# comentário 

plt.show()