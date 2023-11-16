#Vamos a crear una función que me entregue las inversiones de población total y de cada variedad
#para el Hamiltoniano de Jaynes-Cummings y un dipolo de radiación interractuando con un pulso
#Primero veremos una función que muestre la inversión de población para todo el sistema.

import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import math as math

    
#Vamos a crear una función que me entregue las inversiones de población total y de cada variedad
#para el Hamiltoniano de Jaynes-Cummings y un dipolo de radiación interractuando con un pulso
#Primero veremos una función que muestre la inversión de población para todo el sistema.



    



def population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps,tf,tg,retraso):
        #ti va a variar dependiendo del grado de retraso que exista entre heaviside y el pulso
        #El retraso negativo indica que primero actúa g(t) y luego el pulso
        #El retraso negativo minimo indica que pulso y g(t) termianan al mismo tiempo
        #El retraso negativo del 25% indica que cuando termine g(t) al pulso le queda el 25% de su duración
        #El retraso negativo del 50% indica que cuando termine g(t) al pulso le queda el 50% de su duración
        #El retraso negativo del 75% indica que cuando termine g(t) al pulso le queda el 75% de su duración
        #El retraso negativo del 100% indica que cuando termine g(t) al pulso le queda el 100% de su duración

    #Definamos los distintos ti para el inicicio del pulso
    tf1=tf
    if retraso == 0:
        ti = tg-tf
    elif retraso == 25:
        ti = tg-3*tf/4
    elif retraso == 50:
        ti = tg-tf/2
    elif retraso == 75:
        ti = tg-tf/4
    elif retraso == 100:
        ti = tg

    elif retraso == -2:
        ti = tg-2*tf

    elif retraso == -3:
        ti = tg-3*tf
    


    else:
        print("El grado de retraso no es válido")

    #Redefinamos tf
    if retraso == 0:
        tf = ti+tf
    elif retraso == 25:
        tf = ti+tf
    elif retraso == 50:
        tf = ti+tf
    elif retraso == 75:
        tf = ti+tf
    elif retraso == 100:
        tf = ti+tf
    elif retraso == -2:
        tf = ti+tf
    elif retraso == -3:
        tf = ti+tf
    else:
        print("El grado de retraso no es válido")






    t = np.linspace(0,max([tg,tf]),num_steps)
    #t2=tg[-1]

    def heaviside(t,args):
        if t<=tg:
            return np.heaviside(t-0,1)*(g)
        else:
            return 0
        
    
    def heaviside2(t):
        condition = (t > 0) & (t < tg)
        he = np.heaviside(t-0,1)*(g)
        he[~condition] = 0  # Establecer a cero donde la condición no se cumple
        return he
    


    if retraso==0 or retraso==100 or retraso==-2 or retraso==-3:
            
            def pulso(t,args):
                if t>=ti and t<=tf:
                    return E0*(np.sin(np.pi*(t-tf1/4)/(tf-ti))**2) *np.cos(omega_l*(t))
                else:
                    return 0

        
            def pulso2(t):
                condition = (t >= ti) & (t <= tf)
                pulso = E0 * (np.sin(np.pi * (t-tf1/4) / (tf-ti)) ** 2) * np.cos(omega_l * (t-(tf-ti)/2 ))
                pulso[~condition] = 0  # Establecer a cero donde la condición no se cumple
                return pulso


    if retraso==25:
            
            def pulso(t,args):
                if t>=ti and t<=tf:
                    return E0*(np.sin(np.pi*(t-tf1/2)/(tf-ti))**2) *np.cos(omega_l*(t))
                else:
                    return 0

        
            def pulso2(t):
                condition = (t >= ti) & (t <= tf)
                pulso = E0 * (np.sin(np.pi * (t-tf1/2) / (tf-ti)) ** 2) * np.cos(omega_l * (t-(tf-ti)/2 ))
                pulso[~condition] = 0  # Establecer a cero donde la condición no se cumple
                return pulso


    if retraso==50:
            
            def pulso(t,args):
                if t>=ti and t<=tf:
                    return E0*(np.sin(np.pi*(t-3*tf1/4)/(tf-ti))**2) *np.cos(omega_l*(t))
                else:
                    return 0

        
            def pulso2(t):
                condition = (t >= ti) & (t <= tf)
                pulso = E0 * (np.sin(np.pi * (t-3*tf1/4) / (tf-ti)) ** 2) * np.cos(omega_l * (t-(tf-ti)/2 ))
                pulso[~condition] = 0  # Establecer a cero donde la condición no se cumple
                return pulso
            


    if retraso==75:
            
            def pulso(t,args):
                if t>=ti and t<=tf:
                    return E0*(np.sin(np.pi*(t)/(tf-ti))**2) *np.cos(omega_l*(t))
                else:
                    return 0

        
            def pulso2(t):
                condition = (t >= ti) & (t <= tf)
                pulso = E0 * (np.sin(np.pi * (t) / (tf-ti)) ** 2) * np.cos(omega_l * (t-(tf-ti)/2 ))
                pulso[~condition] = 0  # Establecer a cero donde la condición no se cumple
                return pulso

    
    if ini[0]=="e":
        psi0 = tensor(basis(N, ini[1]), basis(2, 1)).unit()

    elif ini[0]=="g":
        psi0 = tensor(basis(N, ini[1]), basis(2, 0)).unit()

    #OPerators
    # operators
    a = tensor(destroy(N), qeye(2))
    sm = tensor(qeye(N), destroy(2))
    #print(a)


    #El Hamiltoniano H0; Cavidad y átomo.
    H0 = omega_c * a.dag() * a + 0.5*omega_0 * (sm.dag() * sm-sm *sm.dag())

    #Hamiltonianos de interacción
    #Primero la parte de interacción del Jaynes-Cummings que tendrá un g(t)
    HI1 = (a.dag()*sm + a*sm.dag())
    #Ahora el dipolo de la radiación
    HI2 = a + a.dag()

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]
    H2=[H0, [HI1,heaviside]]


    # Resolvemos la ecuacion diferencial
    #Incluyamos el operador unidad

    unit= Qobj(tensor(qeye(N), qeye(2)))
    #print(unit)
    

    result1 = mesolve(H, psi0, t, [], [a.dag()*a,(sm.dag() * sm-sm *sm.dag()),sm.dag() * sm,sm *sm.dag(),unit])
    #Necesito el indice para el ultimo tiempo menor o igual a tf
    #Saquemos los estados
    result2 = mesolve(H, psi0, t, [], [])

    inversion=result1.expect[1]
    average_photons=result1.expect[0]
    
    indice = 0
    for i in range(len(t)):
        if t[i]<=ti:
            indice = i
        else:
            break
    
    
    #print(result2.states[indice])
    print(indice)

    #print(result1.expect[0][indice])
    alpha = result1.expect[0][indice]
    coherente=tensor(coherent(N,np.sqrt(alpha)),basis(2,1))
    coherente2=result2.states[indice]
    #print(coherente.norm())
    #print(coherente2.norm())
    #print(coherente.dag()*unit*coherente)
    #print(coherente2.dag()*unit*coherente2)
    print("El número promedio de fotones cuando comienza a actuar g(t) es : ",result1.expect[0][indice])
    #print(result1.expect[4])
    
    #Definamos la función Poissoniana
    x = np.arange(0,N,1)
    def poisson(x, mu):
        return (np.exp(-mu)*mu**x)/math.factorial(int(x))


    result_p=[]

    for i in range(len(x)):
        result_p.append(poisson(x[i],alpha))


    #result_state = mesolve(H, psi0, t, [], [])
        # Resolvemos la ecuacion diferencial en un ciclo for para cada variedad
    #Lista para guardar los resultados
    resultados = []
    resultados_suma = []
    resultadose = []
    resultadosg = []
    """for i in tqdm(range(0,N)):
        if i == 0:
            psig= tensor(basis(N, 0), basis(2, 0)).unit()
            psie= tensor(basis(N, 0), basis(2, 1)).unit()*0
            psie_operator = psie*psie.dag()
            psig_operator = -psig*psig.dag()
            population_operator = psie_operator + psig_operator
            result = mesolve(H, psi0, t, [], [population_operator,psig_operator,psie_operator])
        else:
            psie = tensor(basis(N, i-1), basis(2, 1)).unit()
            psig= tensor(basis(N, i), basis(2, 0)).unit()

            #Definamos sus operadores
            psie_operator = psie*psie.dag()
            psig_operator = -psig*psig.dag()
            population_operator = psie_operator + psig_operator

            #Resolvemos la ecuacion diferencial
            result = mesolve(H, psi0, t, [], [population_operator,psig_operator,psie_operator])
        resultados.append(result.expect[0])
        resultadose.append(result.expect[2])
        resultadosg.append(result.expect[1])
        resultados_suma.append(result.expect[2]-result.expect[1])"""

    b=np.abs(pulso2(t))**2
    unidad=result1.expect[4]
    caja=heaviside2(t)
    pulso_=pulso2(t)
    #Ahora vamos a graficar la inversión de población y el número promedio de fotones
    # Graficamos la inversion de poblacion
    #plt.figure(figsize=(12,7),facecolor='w')
    #plt.grid()
    #plt.title("RWA Jaynes-Cummings y dipolo de la radiacion g={}_E0={} N={}".format(g,E0,N),fontsize=20)
    #plt.plot(t, result1.expect[1], label="Inversion de poblacion",color='#173F5F',lw=2)
    #plt.xticks(np.arange(0, max(t), 200))
    #plt.plot(t, heaviside2(t)*100, label="Envolvente",color='#ED553B',lw=1)
    #plt.plot(t, pulso2(t)*50, label="Pulso",color='green',lw=1)
    #plt.xlabel('Tiempo',fontsize=15)
    #plt.ylabel('Inversion de poblacion',fontsize=15)
    #plt.legend(fontsize=12)
    #plt.savefig('Envolvente_const_retraso/'+"RWA_population_inversion_g={}_E0={}_N={}t={}.png".format(g,E0,N,t[-1]))
    #plt.show()

    # Graficamos el numero promedio de fotones
    #plt.figure(figsize=(12,7))
    #plt.grid()
    #plt.title("RWA Radiación y cavidad g={}_E0={}_N={}{}".format(g,E0,N,area),fontsize=20)
    #plt.plot(t,pulso2(t)*50*30, label="Pulso",color='green',lw=1)
    #plt.plot(t,heaviside2(t)*100*30, label="Envolvente",color='#ED553B',lw=1)
    #plt.plot(t, result1.expect[0], label="Numero promedio de fotones",color='#20639B',lw=2)
    #plt.xlabel('Tiempo',fontsize=15)
    #plt.ylabel('Numero promedio de fotones',fontsize=15)
    #plt.legend(fontsize=12)
    #plt.savefig('Envolvente_const_retraso/'+"RWA_average_photons_g={}_E0={}_N={}_t={}.png".format(g,E0,N,t[-1]))
    #plt.show()



    #Grafiquemos todas las inversiones de población
    """plt.figure(figsize=(12,7))
    plt.grid()
    plt.title("RWA Radiación y cavidad g={}_E0={}_N={}{}".format(g,E0,N,area),fontsize=20)
    for i in range(0,N):
        plt.plot(t, resultados[i], lw=1)

    plt.plot(t, sum(resultados), color="yellow",label="Suma",lw=2)
    plt.plot(t, result1.expect[1], label="Inversion de poblacion",color='#173F5F',lw=1.0)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Inversion de poblacion',fontsize=15)
    plt.legend(fontsize=12)
    #plt.savefig('Envolvente_const_retraso/'+"RWA_population_inversion_all_g={}_E0={}_N={}_t={}.png".format(g,E0,N,t[-1]))
    plt.show()"""



    return t, inversion, average_photons, caja, pulso_
    


   