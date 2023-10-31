#Definamos una poissoiniana
import numpy as np
import matplotlib.pyplot as plt
import math as math


def poisson(x, mu):
    return (np.exp(-mu)*mu**x)/math.factorial(int(x))

#Grafiquemos la función de una partícula en una caja para los primeros 11 niveles
mu = 12.5
x = np.arange(0,50,1)


result=[]

for i in range(len(x)):
    result.append(poisson(x[i],mu))

#Agreguemos más sticks de 10 en 10
x1 = np.arange(0,50,10)
y1 = np.arange(-1,1,0.25)
plt.ylim(-1,1)
plt.xticks(x1)
plt.yticks(y1)

plt.plot(x,result)
plt.show()

