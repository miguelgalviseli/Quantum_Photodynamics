from population_average_const_env_retraso import *
N=40
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=0
n1=1
n2=2
area = "inversion"
ini = ["e", 0]  
num_steps = 5000
tf= (2*n+1)*np.pi/(2*g)+200
tf1= (2*n1+1)*np.pi/(2*g)+200
tf2= (2*n2+1)*np.pi/(2*g)+200
tg=3*tf+1000
retraso=100



#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
print(tg)

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

# population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps, tf, ti, tg)
retraso1 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg, 100)
retraso2 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf1, tg, 100)
retraso3 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf2, tg, 100)




# Graficamos el numero promedio de fotones
plt.figure(figsize=(12,7))
plt.grid()
#plt.title("Inversion de poblacion y numero promedio de fotones",fontsize=20)
plt.plot(retraso1[0], retraso1[2], label="t={}".format(np.round(tf,1)),color='#EF767A')
plt.plot(retraso2[0], retraso2[2], label="t={}".format(np.round(tf1,1)),color='#FFA07A')
plt.plot(retraso3[0], retraso3[2], label="t={}".format(np.round(tf,1)),color='black')

plt.xlabel('Tiempo (u.a.)',fontsize=15)
plt.ylabel(r'$\langle \hat{a}^{\dagger}\hat{a} \rangle$',fontsize=20)
#plt.yticks(np.arange(0, 25, 1))
plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('Aver_plus.png', dpi=250)

plt.show()

#Grafiquemos solamente el pulso y g(t) para retraso = 100%
plt.figure(figsize=(12,7))
plt.grid()

plt.plot(retraso1[0], retraso1[3] * 50, label=r'$\lambda (t)$', color='black', lw=3)
plt.plot(retraso1[0], retraso1[4] * 25, label="pulso", color="#456990", lw=2.5)
plt.xlabel('Tiempo (u.a.)',fontsize=15)
plt.ylabel('Intensidad',fontsize=20)
plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('Pulso_plus.png', dpi=250)

plt.show()