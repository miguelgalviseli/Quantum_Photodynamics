#Vanos a trabajar solamente con la radiación


import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm


def average(N, omega_l,omega_c, E0, ini, num_steps,T):

    t = np.linspace(0,T,num_steps)


    def pulso(t, args):
        return E0*(np.sin(np.pi*t/T)**2) * np.cos(omega_l * (t-T/2))
    
    psi0 = basis(N, ini)

    #OPerators
    a = destroy(N)
    ad = create(N)


    #El Hamiltoniano H0; Cavidad y átomo.
    H0 = omega_c * (ad * a + 0.5) 

    #Hamiltonianos de interacción

    HI1 = (a + ad)

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,pulso]]

    # Resolvemos la ecuacion diferencial
    result = mesolve(H, psi0, t, [], [ad*a])

    #Ahora vamos a graficar la inversión de población y el número promedio de fotones
    # Graficamos la inversion de poblacion
    plt.figure(figsize=(12,7),facecolor='w')
    plt.grid()
    plt.title("Interacción pulso-cavidad",fontsize=20)
    plt.plot(t, result.expect[0], label=r'$\mid \alpha \mid ^{2}$',color='#173F5F',lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Número promedio de fotones',fontsize=15)
    plt.legend(fontsize=15)
    plt.savefig('radiacion_imagenes/' + 'only_radiacion_N={}_E0={}'.format(N,E0) + '.png')
    plt.show()


def average_RWA(N, omega_l,omega_c, E0, ini, num_steps,T):

    t = np.linspace(0,T,num_steps)


    def pulso_plano1(t,args):
        return E0*(np.sin(np.pi*t/T)**2) *np.exp(1j*omega_l*(t-T/2))/2

    #Ahora con argumento negativo
    def pulso_plano2(t,args):
        return E0*(np.sin(np.pi*t/T)**2) *np.exp(-1j*omega_l*(t-T/2))/2
    

    psi0 = basis(N, ini)

    #OPerators
    a = destroy(N)
    ad = create(N)


    #El Hamiltoniano H0; Cavidad y átomo.
    H0 = omega_c * (ad * a + 0.5) 

    #Hamiltonianos de interacción

    HI1 = a
    HI2 = ad 
    H = [H0, [HI1, pulso_plano1], [HI2, pulso_plano2]]


    # Resolvemos la ecuacion diferencial
    result = mesolve(H, psi0, t, [], [ad*a])

    #Ahora vamos a graficar la inversión de población y el número promedio de fotones
    # Graficamos la inversion de poblacion
    plt.figure(figsize=(12,7),facecolor='w')
    plt.grid()
    plt.title("Interacción pulso-cavidad con RWA",fontsize=20)
    plt.plot(t, result.expect[0], label=r'$\mid \alpha \mid ^{2}$',color='#173F5F',lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Número promedio de fotones',fontsize=15)
    plt.legend(fontsize=15)
    plt.savefig('radiacion_imagenes/' + 'only_radiacion_N={}_E0={}_RWA'.format(N,E0) + '.png')
    plt.show()
