from population_average_const_env_retraso_negative import *
N=40
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
#n=2
#Voy  dejar fijo n para que el tiempo de pulso siempre sea el mismo
n=2
n1=6
n2=3
n3=6
tf= (2*n+1)*np.pi/(2*g)
area = "inversion"
ini = ["e", 0]  
num_steps = 5000
#tf= (2*n+1)*np.pi/(2*g)
#tg=3*tf+1000
tg1=(n1)*np.pi/(2*g)
tg2=(2*n2+1)*np.pi/(2*g)
tg3=(2*n3+1)*np.pi/(4*g)

retraso=100



#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
#print(tg)

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

# population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps, tf, ti, tg)
retraso1 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg1, 100)
retraso2 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg2, 100)
retraso3 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg3, 100)


# Crear un subplot con 5 filas y 1 columna
fig, axs = plt.subplots(3, 1, figsize=(16, 18), facecolor='w', sharex=True)
# Ajustar el tamaño de la figura
fig.subplots_adjust(hspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.05)

# Graficar cada retrazo
for i, retraso in enumerate([retraso1, retraso2, retraso3]):
    #En cada esquina superior izquierda de cada subplot, agregar A) B) C) D) E)
    axs[i].text(0, 0.90,  chr(65 + i) + ")", transform=axs[i].transAxes, size=20, weight='bold')
    axs[i].plot(retraso[0], retraso[1], label=r'$\langle \hat{\sigma}_{3} \rangle$', color="#EF767A", lw=2.5)
    axs[i].plot(retraso[0], retraso[3] * 50, label=r'$\lambda (t)$', color='black', lw=3)
    axs[i].plot(retraso[0], retraso[4] * 25, label="pulso", color="#C70039", lw=2)
    axs[i].tick_params(axis='y', labelsize=18)
    axs[i].legend(fontsize=18)
    axs[i].grid()

    # Etiquetas en x solo para el último subplot
    if i == 2:
        axs[i].set_xlabel("Tiempo (u.a)", fontsize=18)
        axs[i].set_xticks([valor for valor in axs[i].get_xticks() if valor >= axs[i].get_xlim()[0] and valor <= axs[i].get_xlim()[1]])
        axs[i].tick_params(axis='x', labelsize=18)
        axs[i].tick_params(axis='y', labelsize=18)

    # Etiquetas en x solo para el último subplot
    if i == 0 or i == 1 or i == 2:
         axs[i].plot(retraso[0], retraso[1], label=r'$\langle \hat{\sigma}_{3} \rangle$ Analítico', color="blue", lw=0.75)
         axs[i].legend(fontsize=18)

# Unificar leyendas
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.1, 0.90), fontsize=15)

# Ajustar el diseño
plt.tight_layout()

# Ajustar el tamaño de la imagen al guardarla (aumentando el dpi)
plt.savefig('Inversiones_minus_3.png', dpi=250)

# Mostrar el gráfico
plt.show()

# Graficamos el numero promedio de fotones
plt.figure(figsize=(12,7))
plt.grid()
#plt.title("Inversion de poblacion y numero promedio de fotones",fontsize=20)
plt.plot(retraso1[0], retraso1[2], label="Terminando en |e,0>",color='#3CAEA3')
plt.plot(retraso2[0], retraso2[2], label="Terminando en |g,1>",color='black')
plt.plot(retraso3[0], retraso3[2], label=r'$\frac{1}{\sqrt{2}}(|e,0\rangle + |g,1\rangle)$',color='#ED553B')
#Quiero poner más stickers en el eje y, desde 0 hsta 18, de 1 en 1
plt.yticks(np.arange(0, 18, 1))
plt.xlabel('Tiempo (u.a.)',fontsize=15)
plt.ylabel(r'$\langle \hat{a}^{\dagger}\hat{a} \rangle$',fontsize=15)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#Quiero poner una linea punteada en los valores de cada valor del número promedio de fotones
#El valor de y es el último valor de retraso1[2], retraso2[2] y retraso3[2]
plt.hlines(retraso1[2][-1], retraso1[0][-1], 2000, linestyles='dashed', colors='black')
plt.hlines(retraso2[2][-1], retraso2[0][-1], 2000, linestyles='dashed', colors='black')
plt.hlines(retraso3[2][-1], retraso3[0][-1], 2000, linestyles='dashed', colors='black')
#Quiero un solo legend para las lineas punteadas
plt.plot([], [], 'k--', label='Valor final para cada caso', lw=2)
plt.legend(fontsize=16)
plt.savefig('Aver_minus-3.png', dpi=250)

plt.show()
