import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import cv2
import os
from population_average_const_env_retraso_negative import *

import math as math
N=50
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
area = "inversion"
ini = ["e", 0]  
num_steps = 1000

#Tiempos retraso máximo
n=2
tg= (2*n+1)*np.pi/(2*g)
ti=tg
tf=2*tg
#Tiempos retraso cero
#n=3
#tf= (2*n+1)*np.pi/(2*g)
#tg=tf
#ti=0


#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print("El tiempo de retraso es {}".format(round(tf-tg,2)))
print("El tiempo de pulso es {}".format(round(tf-tg,2)))
print("El tiempo de g(t) es {}".format(round(tg,2)))
print("El tiempo total es {}".format(round(ti+tf,2)))

#print(T)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
t = np.linspace(0,max([tf,tg+ti]),num_steps)
#Tiempo para los valores teoricos en el mismo intervalo que los valores simulados


a = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg)
t2  = np.linspace(0,tg,len(t)-a[6]-1)
 #Definamos las funciones teóricas de la inversión de población oscilando

def poisson_e(mu,n,t,g):
    return ((np.exp(-mu)*mu**n)/math.factorial(n))*(np.cos(g*np.sqrt(n+1)*t))**2

def poisson_g(mu,n,t,g):    
    return -((np.exp(-mu)*mu**n)/math.factorial(n))*(np.sin(g*np.sqrt(n+1)*t))**2

#Ahora guardemos los tiempos despues del indice
t_sim=[]
for i in range(len(a[0][1])):
    if i >a[6]:
        t_sim.append(t[i])

#Volvamos t_sim un array
t_sim=np.array(t_sim)

#Necesitamos guardar los valores a partir de un indice
e0_theoric=[]
g1_theoric=[]
e0=poisson_e(a[5],18,t2,g)
g1=poisson_g(a[5],18,t2,g)

print(a[6])





#Ahora e0 y g1 simulados despues del indice
e0_sim=[]
g1_sim=[]
for i in range(len(a[0][1])):
    if i >a[6]:
        e0_sim.append(a[0][19][i])
        g1_sim.append(a[1][19][i])





#Grafiquemos la función teórica para la primera variedad de la inversión de población y los valores de la simulación

plt.figure(figsize=(12,7))
plt.plot(t_sim,e0, color="#445D48",label="En |e,0> teórico",linewidth=3)
plt.plot(t_sim,g1,color="#5E3023",label="En |g,1> teórico",linewidth=3)
plt.plot(t_sim,e0_sim,color="#402149",label="En |e,0> simulado",linewidth=0.9)
plt.plot(t_sim,g1_sim,color="#A41F13",label="En |g,1> simulado",linewidth=0.9)
plt.xticks(np.arange(0, max(t), 200))
plt.grid()
plt.xlabel("Tiempo")
plt.ylabel("Poblaciones")
plt.legend()
plt.show()

# Graficamos la unidad
plt.figure(figsize=(12,7))

plt.title("Operador identidad",fontsize=20)

plt.plot(t, a[8], label="Unidad",color='red')
plt.grid()
plt.xlabel('Tiempo',fontsize=15)
plt.ylabel('Unidad',fontsize=15)
plt.ylim(0,1.5)
plt.legend(fontsize=12)
plt.show()



#print(t_sim)
