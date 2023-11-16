from intermedio_population_average_const_env_retraso import *
N=40
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=2
area = "inversion"
ini = ["e", 0]  
num_steps = 5000
tf= (2*n+1)*np.pi/(2*g)
tg=3*tf+1000
retraso=100



#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
print(tg)

import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

# population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps, tf, ti, tg)
retraso5 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg, -2)
retraso4 = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, area, ini, num_steps, tf, tg, -3)


# Crear un subplot con 5 filas y 1 columna
fig, axs = plt.subplots(2, 1, figsize=(12, 14), facecolor='w', sharex=True)
# Ajustar el tamaño de la figura
fig.subplots_adjust(hspace=0.5, left=0.1, right=0.9, top=0.95, bottom=0.05)

# Graficar cada retrazo
for i, retraso in enumerate([ retraso4, retraso5]):
    #En cada esquina superior izquierda de cada subplot, agregar A) B) C) D) E)
    axs[i].text(0, 0.90,  chr(65 + i) + ")", transform=axs[i].transAxes, size=20, weight='bold')
    axs[i].plot(retraso[0], retraso[1], label=r'$\langle \hat{\sigma}_{3} \rangle$', color="#FFC30F", lw=2.5)
    axs[i].plot(retraso[0], retraso[3] * 50, label=r'$\lambda (t)$', color='black', lw=3)
    axs[i].plot(retraso[0], retraso[4] * 25, label="pulso", color="#C70039", lw=2)
    axs[i].tick_params(axis='y', labelsize=18)
    axs[i].legend(fontsize=18)
    axs[i].grid()

    # Etiquetas en x solo para el último subplot
    if i == 1:
        axs[i].set_xlabel("Tiempo (u.a)", fontsize=18)
        axs[i].set_xticks([valor for valor in axs[i].get_xticks() if valor >= axs[i].get_xlim()[0] and valor <= axs[i].get_xlim()[1]])
        axs[i].tick_params(axis='x', labelsize=18)
        axs[i].tick_params(axis='y', labelsize=18)



# Unificar leyendas
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.1, 0.90), fontsize=15)

# Ajustar el diseño
plt.tight_layout()

# Ajustar el tamaño de la imagen al guardarla (aumentando el dpi)
plt.savefig('Inversiones_minus_int.png', dpi=200)

# Mostrar el gráfico
plt.show()

# Graficamos el numero promedio de fotones
plt.figure(figsize=(12,7))
plt.grid()
#plt.title("Inversion de poblacion y numero promedio de fotones",fontsize=20)
plt.plot(retraso5[0], retraso5[2], label="Retraso = -300%",color='#173F5F')
plt.plot(retraso4[0], retraso4[2], label="Retraso = -200%",color='#FFC30F')

plt.xlabel('Tiempo (u.a.)',fontsize=15)
plt.ylabel(r'$\langle \hat{a}^{\dagger}\hat{a} \rangle$',fontsize=15)
plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('Aver_minus_int.png', dpi=250)

plt.show()