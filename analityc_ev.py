import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import math as math
from scipy.special import genlaguerre
from scipy.special import factorial
from population_average_const_env_retraso_negative import *

#Veamos la inversión de población para un pulso actuando sobre un estado inicial venido del Jaynes-Cummings
N=50
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=2
area = "inversion"
ini = ["e", 0]  
num_steps = 1000
tg= (n)*np.pi/(g)
tf=2*tg
ti=tg
t = np.linspace(0,max([tf,tg+ti]),num_steps)

#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
print(ti)
print(tg)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg)
a=population_inversion_all(N, omega_l, omega_0, omega_c,g, E0, n, area, ini, num_steps,tf,ti,tg)
print(a[10])
print(a[9])
#Definamos las funciones que calculan los elementos de matriz en la expresión analítica

#Elementos de matriz para l mayor o igual que n
def matrix_element(m,n, alpha):
        if m>=n:
            n=float(n)
            m=float(m)
            alpha=float(alpha)
            return (alpha**(m-n))*np.sqrt(factorial(n)/factorial(m))*np.exp(-(abs(alpha)**2)/2)*genlaguerre(n,m-n)(abs(alpha)**2)
        else:
             return ((-alpha)**(n-m))*np.sqrt(factorial(m)/factorial(n))*np.exp(-(abs(alpha)**2)/2)*genlaguerre(m,n-m)(abs(alpha)**2)
        
             

#Definamos la función que calcula la inversión de población analítica

def inversion_analitica(t, N, alpha,g):
    list=[]
    for k in tqdm(range(len(t))):
        suma=0
        for i in range(N):
            for j in range(N):
                suma+= abs(matrix_element(i,j,np.sqrt(alpha[k])))**2*np.cos(g*np.sqrt(i+1)*t[k])**2 - abs(matrix_element(i,j+1,np.sqrt(alpha[k])))**2*np.sin(g*np.sqrt(i+1)*t[k])**2
        list.append(suma/20)
    return list


t_sim=[]
for i in range(len(a[0][1])):
    if i >a[6]:
        t_sim.append(t[i])

#Volvamos t_sim un array
t_sim=np.array(t_sim)
#Tiempo para la ecuación
t2  = np.linspace(0,tg,len(t)-a[6]-1)

#print(t_sim)
#print(t2)
#Extraigamos los alphas que me interesan
alpha=[]
for i in range(len(a[0][1])):
    if i >a[6]:
        alpha.append(a[10][i])
print(len(alpha))

#eXTRAIGAMOS LAS INVERSIONES DE POBLACIÓN QUE ME INTERESAN
inversion=[]
for i in range(len(a[0][1])):
    if i >a[6]:
        inversion.append(a[9][i])

#GRAFICAMOS
c=inversion_analitica(t2,N,alpha,g)
#print(c)
plt.plot(t_sim,inversion, label="Simulación")
plt.plot(t_sim,c, label="Teórico")
plt.legend()
plt.show()

                  