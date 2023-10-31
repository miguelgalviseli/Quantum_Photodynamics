#Vamos a crear una función que me entregue las inversiones de población total y de cada variedad
#para el Hamiltoniano de Jaynes-Cummings y un dipolo de radiación interractuando con un pulso
#Primero veremos una función que muestre la inversión de población para todo el sistema.

import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm

    

    



def population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg):


    t = np.linspace(0,max([tg,tf]),num_steps)

    def heaviside(t,args):
        return (np.heaviside(t-0,1)-np.heaviside(t-tg,1))*g
    
    def heaviside2(t):
        condition = (t >= 0) & (t <= tg)
        heaviside = (np.heaviside(t-0,1)-np.heaviside(t-tg,1))*g
        heaviside[~condition] = 0
        return heaviside

    def pulso(t,args):
        if t>=ti and t<=tf:
            return E0*(np.sin(np.pi*t/(tf-ti))**2) *np.cos(omega_l*(t-(tf-ti)/2))
        else:
            return 0
        
    def pulso2(t):
        condition = (t >= ti) & (t <= tf)
        pulso = E0 * (np.sin(np.pi * t / (tf-ti)) ** 2) * np.cos(omega_l * (t - (tf-ti) / 2)) 
        pulso[~condition] = 0  # Establecer a cero donde la condición no se cumple
        return pulso




    
    if ini[0]=="e":
        psi0 = tensor(basis(N, ini[1]), basis(2, 1))

    elif ini[0]=="g":
        psi0 = tensor(basis(N, ini[1]), basis(2, 0))

    #OPerators
    # operators
    a = tensor(destroy(N), qeye(2))
    sm = tensor(qeye(N), destroy(2))


    #El Hamiltoniano H0; Cavidad y átomo.
    H0 = omega_c * a.dag() * a + 0.5*omega_0 * (sm.dag() * sm-sm *sm.dag())

    #Hamiltonianos de interacción
    #Primero la parte de interacción del Jaynes-Cummings que tendrá un g(t)
    HI1 = (a.dag()*sm + a*sm.dag())
    #Ahora el dipolo de la radiación
    HI2 = a + a.dag()

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]


    # Resolvemos la ecuacion diferencial
    result1 = mesolve(H, psi0, t, [], [a.dag()*a,(sm.dag() * sm-sm *sm.dag()),sm.dag() * sm,sm *sm.dag()])
    #result_state = mesolve(H, psi0, t, [], [])
        # Resolvemos la ecuacion diferencial en un ciclo for para cada variedad
    #Lista para guardar los resultados
    resultados = []
    resultadose = []
    resultadosg = []
    for i in tqdm(range(0,N)):
        if i == 0:
            psig= tensor(basis(N, 0), basis(2, 0))
            psig_operator = -psig*psig.dag()
            result = mesolve(H, psi0, t, [], [psig_operator,psig_operator,psig_operator])
        else:
            psie = tensor(basis(N, i-1), basis(2, 1))
            psig= tensor(basis(N, i), basis(2, 0))

            #Definamos sus operadores
            psie_operator = psie*psie.dag()
            psig_operator = -psig*psig.dag()
            population_operator = psie_operator + psig_operator

            #Resolvemos la ecuacion diferencial
            result = mesolve(H, psi0, t, [], [population_operator,psig_operator,psie_operator])
        resultados.append(result.expect[0])
        resultadose.append(result.expect[2])
        resultadosg.append(result.expect[1])
    b=np.abs(pulso2(t))**2
    #Ahora vamos a graficar la inversión de población y el número promedio de fotones
    # Graficamos la inversion de poblacion
    plt.figure(figsize=(12,7),facecolor='w')
    plt.grid()
    plt.title("RWA Jaynes-Cummings y dipolo de la radiacion g={}_E0={} N={}".format(g,E0,N),fontsize=20)
    plt.plot(t, result1.expect[1], label="Inversion de poblacion",color='#173F5F',lw=2)
    plt.plot(t, heaviside2(t)*100, label="Envolvente",color='#ED553B',lw=1)
    plt.plot(t, pulso2(t)*50, label="Pulso",color='green',lw=1)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Inversion de poblacion',fontsize=15)
    plt.legend(fontsize=12)
    plt.savefig('Envolvente_const_retraso/'+"RWA_population_inversion_g={}_E0={}_N={}t={}.png".format(g,E0,N,t[-1]))
    plt.show()

    # Graficamos el numero promedio de fotones
    plt.figure(figsize=(12,7))
    plt.grid()
    plt.title("RWA Radiación y cavidad g={}_E0={}_N={}{}".format(g,E0,N,area),fontsize=20)
    plt.plot(t,pulso2(t)*50*30, label="Pulso",color='green',lw=1)
    plt.plot(t,heaviside2(t)*100*30, label="Envolvente",color='#ED553B',lw=1)
    plt.plot(t, result1.expect[0], label="Numero promedio de fotones",color='#20639B',lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Numero promedio de fotones',fontsize=15)
    plt.legend(fontsize=12)
    plt.savefig('Envolvente_const_retraso/'+"RWA_average_photons_g={}_E0={}_N={}_t={}.png".format(g,E0,N,t[-1]))
    plt.show()

    return resultadose, resultadosg, b
    


   