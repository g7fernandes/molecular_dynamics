 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def unprecize(step,p1,p2,sigma,epsil):
    r = np.linspace(p1,p2,int((p2-p1)/step))
    f = (1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil*r
    I = 0
    for i in range(len(f)-1):
        I += 0.5*(f[i+1]+f[i])*step
    return I

def trapz(x,y):
    I = 0
    for i in range(len(x)-1):
        I += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i-1])
    return I


def flj(r,sigma, epsil):
    return (1/r**2)*(sigma/r)**6*(1-2*(sigma/r)**6)*24*epsil*r

print("O que fazer?\n1 - Integrar LJ trapezoidal.\n2 - Distribuição Maxwell-Boltzmann\n3 - Interagir duas partículas.\n")
tarefa = input("Entre 1, 2 ou 3: ")

if tarefa == "1":

    v_max = 5 
    n_points = 100
    dt = np.logspace(-5,-1,n_points)
    step = v_max*dt
    sigma = 1
    epsil = 1
    m = 1
    p1,p2 = sigma*.98, 3*sigma
    desv = np.zeros(n_points)

    r = np.linspace(p1,p2,10000)
    f = flj(r,sigma,epsil)
    v = 4*epsil*((sigma/r)**6*((sigma/r)**6)-1)

    E_precisa = quad(flj,p1,p2,args=(sigma,epsil))

    for i in range(n_points):
        desv[i] = (unprecize(step[i],p1,p2,sigma,epsil)-E_precisa[0])/E_precisa[0]


    fig, ax1 = plt.subplots()
    ax2  = ax1.twinx()

    dts = [0.00001, 0.0001, 0.001, 0.01, 0.02]

    lns2 = ax2.plot(r,0.00001*f/m,label='v after dt = 0.00001')
    lns3 = ax2.plot(r,0.0001*f/m,label='v after dt = 0.0001')
    lns4 = ax2.plot(r,0.001*f/m,label='v after dt = 0.001')
    lns5 = ax2.plot(r,0.01*f/m,label='v after dt = 0.01')
    lns6 = ax2.plot(r,0.01*f/m,label='v after dt = 0.02')
    #lns7 = ax2.plot(r,0.05*f/m,label='v after dt = 0.05')
    ax2.set_ylabel('speed')

    lns1 = ax1.plot(r,f,'-k',label='force')
    ax1.set_xlabel('r')
    ax1.set_ylabel('force')

    leg = lns2 + lns3 + lns4 + lns5 + lns6 #+ lns7
    labs = [l.get_label() for l in leg]
    ax1.legend(leg, labs, loc=0)

    plt.show()

    plt.semilogx(dt,desv,label='desvio')
    plt.legend()

    plt.show()
elif tarefa == "2":

# pdf maxwell

    plt.figure(7)
    ans = "y"
    while (ans == "y"):

        prop = input("Enter Temperature and mass and v_max for the graph: ").split()
        T = float(prop[0])
        m = float(prop[1])

        print("mean v = {} at T + {}".format(np.sqrt(2*T),T))

        v = np.linspace(0,float(prop[2]),1000)

        pdf = (m/(2*np.pi*T))**(3/2)*(4*np.pi*v**2)*np.exp(-(m*v**2/2*T))

        lab = "T = " + prop[0] + " m = " + prop[1]
        plt.plot(v,pdf,label=lab)
        plt.xlabel("Velocity")
        plt.ylabel("Probabilty density")
        ans = input("Draw one more curve? [y/n]")
    plt.legend()
    plt.show()

elif tarefa == 3:

# caso para duas partículas muito rápidas se aproximando


    t_fim  = 2
    print("tfim = {}\n".format(t_fim))
    vini = -float(input('V inicial >0 ? '))

    p1 =  3*sigma/2
    p2 =  -3*sigma/2
    ch = 'c'
    w = np.zeros(len(dts))
    j = 0
    while ch != 'cq':
        for dt1 in dts:

            v1, v2 = np.zeros(int(t_fim/dt1)), np.zeros(int(t_fim/dt1))
            x1, x2 = np.zeros(int(t_fim/dt1)), np.zeros(int(t_fim/dt1))
            F = np.zeros(int(t_fim/dt1))
            x1[0] = p1
            x2[0] = p2
            v1[0] = vini
            v2[0] = -vini
            for i in range(len(v1)-1):
                f1 = -flj(x1[i]-x2[i],sigma,epsil)

                x1[i+1] = x1[i] + dt1*v1[i] + f1*dt1**2/(2*m)
                x2[i+1] = x2[i] + dt1*v2[i] + (-f1)*dt1**2/(2*m)
                
                f2 = -flj(x1[i+1]-x2[i+1],sigma,epsil)
                F[i] = f2

                v1[i+1] = v1[i] + (f1+f2)*dt1/(2*m) 
                v2[i+1] = v2[i] -(f1+f2)*dt1/(2*m)

            temp = np.linspace(0,t_fim, int(t_fim/dt1))
            plt.figure(3)
            plt.plot(temp,v1,label='V1 dt = '+str(dt1))
            #plt.plot(temp,v2,label='V2 dt = '+str(dt1))
            plt.figure(4)
            plt.plot(temp,x1,label='X1 dt = '+str(dt1))
            plt.figure(5)
            plt.plot(temp,F,label='F1 dt '+str(dt1))
            w[j] = trapz(x1,F)
            j += 1

        plt.figure(4)
        plt.xlabel('Time')
        plt.ylabel('Distance')
        plt.legend()

        plt.figure(3)
        plt.xlabel('Time')
        plt.ylabel('Speed')
        plt.legend()

        plt.figure(5)
        plt.xlabel('temp')
        plt.ylabel('force')
        plt.legend()

        plt.figure(6)
        plt.semilogx(dts,w)
        plt.xlabel('dt')
        plt.xlabel('work on one particle')
        
        plt.show()

        ch = input("Change something? \nEnter vini for initial velocity;\nEnter p1 for initial position,\nEnter tfim for final time.\n> Enter c to continue or cq to quit.\n")
        while (ch[0] != 'c'):
            if ch == 'vini':
                vini = -float(input('-V inicial?  '))
            elif ch == 'p1':
                vini = float(input('p1 inicial? '))
            elif ch == 'tfim':
                t_fim = float(input('tfim ? '))
            ch = input("Change something? \nEnter vini for initial velocity;\nEnter p1 for initial position,\nEnter tfim for final time.\n> Enter c to continue or cq to quit.\n")

            

    #plt.plot(temp,x1,label='X1 dt = '+str(dt1))
    #plt.plot(temp,x2,label='X2 dt = '+str(dt1))



