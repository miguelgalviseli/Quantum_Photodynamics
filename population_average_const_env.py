#Vamos a crear una función que me entregue las inversiones de población total y de cada variedad
#para el Hamiltoniano de Jaynes-Cummings y un dipolo de radiación interractuando con un pulso
#Primero veremos una función que muestre la inversión de población para todo el sistema.

import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm


def population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps):
    if area=="inversion":
        T = (2*n+1)*np.pi/(g)

    elif area=="recuperacion":
        T = np.pi*n/(g)

    elif area=="zero":
        T=(2*n+1)*np.pi/(2*g)

    t = np.linspace(0,T,num_steps)

    def heaviside(t,args):
        return g
    
    def heaviside2(t,g):
        return g
    def pulso(t, args):
        return E0*(np.sin(np.pi*t/T)**2) * np.cos(omega_l * (t-T/2))
    
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
    HI2 = (a + a.dag())

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]

    # Resolvemos la ecuacion diferencial
    result = mesolve(H, psi0, t, [], [a.dag()*a,(sm.dag() * sm-sm *sm.dag()),sm.dag() * sm,sm *sm.dag()])

    #Ahora vamos a graficar la inversión de población y el número promedio de fotones
    # Graficamos la inversion de poblacion
    plt.figure(figsize=(12,7),facecolor='w')
    plt.grid()
    plt.title("Jaynes-Cummings y dipolo de la radiacion g={} E0={} N={}, {}".format(g,E0,N,area),fontsize=20)
    plt.plot(t, result.expect[1], label="Inversion de poblacion",color='#173F5F',lw=2)
    plt.plot(t, result.expect[2], label="|e>",color='#F6D55C',lw=2)
    plt.plot(t, result.expect[3], label="|g>",color='#ED5538',lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Inversion de poblacion',fontsize=15)
    plt.title("Jaynes-Cummings y pulso con RWA",fontsize=20)
    plt.legend(fontsize=15)
    plt.savefig("population_inversion_general/population_inversion_g={}_E0={}_N={} {}.png".format(g,E0,N,area))
    plt.show()
    

    # Graficamos el numero promedio de fotones
    plt.figure(figsize=(12,7))
    plt.grid()
    plt.title("Jaynes-Cummings y dipolo de la radiacion g={} E0={} N={}".format(g,E0,N),fontsize=20)
    plt.plot(t, result.expect[0], label="Numero promedio de fotones",lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Numero promedio de fotones',fontsize=15)
    plt.title("Jaynes-Cummings y pulso con RWA",fontsize=20)
    plt.legend(fontsize=15)
    plt.savefig("average_photons_general/average_photons_g={}_E0={}_N={} {}.png".format(g,E0,N,area))
    plt.show()
    


        
#Ahora veremos una función que muestre la inversión de población para cada variedad pedida 


def population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps):
    if area=="inversion":
        T = (2*n+1)*np.pi/(g)

    elif area=="recuperacion":
        T = np.pi*n/g

    elif area=="zero":
        T=(2*n+1)*np.pi/(2*g)

    t = np.linspace(0,T,num_steps)

    def heaviside(t,args):
        return g
    
    def heaviside2(t,g):
        return g
    def pulso(t, args):
        return E0*(np.sin(np.pi*t/T)**2) * np.cos(omega_l * (t-T/2))
    
    def pulso_(t):
        return abs(0.2*np.cos(omega_l * (t-T/2)))**2
    
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
    HI2 = (a + a.dag())

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]

    # Resolvemos la ecuacion diferencial en un ciclo for para cada variedad
    #Lista para guardar los resultados
    resultados = []
    resultadose = []
    resultadosg = []
    for i in tqdm(variedades):
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

    #Linea prueba para constatar que la suma está bien hecha cuando se tienen todas las variedades
    #result2 = mesolve(H, psi0, t, [], [a.dag()*a,(sm.dag() * sm-sm *sm.dag()),sm.dag() * sm,sm *sm.dag()])
    #Ahora vamos a graficar la inversión de población de las variedades pedidas
    # Graficamos la inversion de poblacion
    plt.figure(figsize=(30,16),facecolor='w')
    plt.suptitle("Jaynes-Cummings y dipolo de la radiacion",fontsize=20)

    
    for i in range(len(variedades)):
        plt.subplot(4,3,i+1)
        plt.bar(t, resultadose[i], label='|e>'+str(variedades[i]),lw=2)
        plt.bar(t, resultadosg[i], label='|g>'+str(variedades[i]),lw=2)
        plt.legend(fontsize=10)
        #plt.plot(t,pulso_(t),label="Pulso",color="red",lw=2)
    #plt.plot(t, sum(resultados), label="Suma de las inversiones",lw=2)
    #plt.plot(t, result2.expect[1], label="Inversion de poblacion general",color='#173F5F',lw=2)
    #plt.xlabel('Tiempo',fontsize=15)
    #plt.ylabel('Inversion de poblacion',fontsize=15)
    #plt.title("Jaynes-Cummings y pulso con RWA",fontsize=20)
    plt.savefig("Population_inversion_variedades.png")
    plt.show()
 
    


            
            
#Ahora una función para el número promedio de fotones por variedad

def average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps):
    if area=="inversion":
        T = (2*n+1)*np.pi/(g)

    elif area=="recuperacion":
        T = np.pi*n/g

    elif area=="zero":
        T=(2*n+1)*np.pi/(2*g)

    t = np.linspace(0,T,num_steps)

    def heaviside(t,args):
        return g
    
    def heaviside2(t,g):
        return g
    def pulso(t, args):
        return E0*(np.sin(np.pi*t/T)**2) * np.cos(omega_l * (t-T/2))
    
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
    HI2 = (a + a.dag())

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]

    # Resolvemos la ecuacion diferencial
    N_=a.dag()*a
    N_matriz=N_.full()
    Lista_expect=[]  #Lista para guardar los resultados    
    probabilidad=[] #Lista para guardar las probabilidades de cada variedad
    #Quiero hacer un histograma para cada tiempo en n
    for i in tqdm(range(0,N)):
        matriz_1=np.zeros((2*N,2*N)) #Matriz para el operador número de fotones en la variedad i
        matriz_1[2*i-1][2*i-1] =1
        matriz_1[2*i][2*i] =1
        matriz = [[0] * N*2 for _ in range(N*2)]
        matriz[2*i-1][2*i-1] = 1
        matriz[2*i][2*i] = 1
        producto=matriz*N_matriz
        p_i=Qobj(matriz_1,dims=[[N,2],[N,2]])
        producto=Qobj(producto,dims=[[N,2],[N,2]])
        result1 = mesolve(H, psi0, t,[], [producto,p_i])
        Lista_expect.append(result1.expect[0])
        probabilidad.append(result1.expect[1])
    


    #Grafiquemos la suma y cada uno de los elementos de la lista
    result2 = mesolve(H, psi0, t, [], [a.dag()*a,(sm.dag() * sm-sm *sm.dag()),sm.dag() * sm,sm *sm.dag()])
    plt.figure(figsize=(12,7),facecolor='w')
    plt.grid()
    plt.title("Jaynes-Cummings y dipolo de la radiacion",fontsize=20)
    for i in range(0,N):
        plt.plot(t, Lista_expect[i], label="Variedad "+str(i),lw=2)
    plt.plot(t, sum(probabilidad), label="Suma de los fotones",lw=1.5)
    plt.plot(t, result2.expect[0], label="Numero promedio de fotones",color='red',lw=2)
    plt.xlabel('Tiempo',fontsize=15)
    plt.ylabel('Numero promedio de fotones',fontsize=12)
    #plt.legend(fontsize=10)
    plt.savefig("Average_photons_individual.png")
    plt.show()







    
    return probabilidad
    

    
        

    

#Ahora veremos una función que muestre la inversión de población para cada variedades 


def population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps):
    if area=="inversion":
        T = (2*n+1)*np.pi/(g)

    elif area=="recuperacion":
        T = np.pi*n/g

    elif area=="zero":
        T=(2*n+1)*np.pi/(2*g)

    t = np.linspace(0,T,num_steps)

    def heaviside(t,args):
        return g
    
    def heaviside2(t,g):
        return g
    def pulso(t, args):
        return E0*(np.sin(np.pi*t/T)**2) * np.cos(omega_l * (t-T/2))
    

    
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
    HI2 = (a + a.dag())

    #Definimos el Hamiltoniano total que será usado por Qutip

    H = [H0, [HI1,heaviside],[HI2,pulso]]

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
    return resultadose, resultadosg

