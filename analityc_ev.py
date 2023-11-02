import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import math as math
from scipy.special import genlaguerre
from scipy.special import factorial
from population_average_const_env_retraso_negative import *

#Veamos la inversión de población para un pulso actuando sobre un estado inicial venido del Jaynes-Cummings
N=20
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=5
area = "inversion"
ini = ["e", 0]  
num_steps = 1000
tg= (n)*np.pi/(2*g)
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
        list.append(suma)
    return list

    

A=inversion_analitica(t, N, a[10],g)

plt.plot(t,A)
plt.plot(t,a[9])
plt.show()
                  