#Definamos la funcion de una partícula en una caja
import numpy as np
import matplotlib.pyplot as plt

def box(x, L,n):
    return np.sqrt(2/L)*np.sin(n*np.pi*x/L)

#Grafiquemos la función de una partícula en una caja para los primeros 11 niveles
L = 11
x = np.linspace(0,L,1000)

for i in range(1,11):
    plt.plot(x,box(x,L,i)+i, label = "n = "+str(i))
    plt.hlines(i,0,L,linestyles = "dashed")
plt.legend()
plt.show()





