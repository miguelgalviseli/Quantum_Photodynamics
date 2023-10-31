from population_average_const_env_retraso import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import os

N=75
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=10
area = "inversion"
ini = ["e", 0]  
num_steps = 500
tf= (2*n+1)*np.pi/(2*g)-1500
tg=3*tf
ti=tf

#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
print(ti)
print(tg)


t = np.linspace(0,max([tf,tg+ti]),num_steps)

a = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg)

# Directorio para guardar los archivos temporales
output_directory = 'temp_images'

# Crear el directorio si no existe
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Función para guardar los histogramas en archivos temporales
def save_histogram(e, g, tiempo, t):
    lista_enelementose = [item[tiempo] for item in e]
    lista_enelementosg = [item[tiempo] for item in g]
    lista_n = list(range(len(e)))

    plt.figure(figsize=(10, 6))
    plt.bar(lista_n, lista_enelementose, color='#173F5F', edgecolor='black', width=1, label='Estate')
    plt.bar(lista_n, lista_enelementosg, color='orange', edgecolor='black', width=1, label='Ground State')
    plt.xlabel('Variedad n', fontsize=15)
    plt.ylabel('Inversión', fontsize=15)
    plt.title(f'Jaynes-Cummings y radiación. t = {np.round(t[tiempo], 2)}', fontsize=15)
    plt.axhline(y=1, color='g', linestyle='-')
    plt.axhline(y=-1, color='g', linestyle='-')
    plt.legend(loc='upper right')
    
    # Guardar la figura en un archivo temporal
    image_path = os.path.join(output_directory, f'hist_{tiempo}.png')
    plt.savefig(image_path)
    plt.close()

# Crear una lista de nombres de archivos de imágenes
image_files = [os.path.join(output_directory, f'hist_{i}.png') for i in range(num_steps)]

# Función para eliminar archivos temporales
def cleanup_files(files):
    for file in files:
        os.remove(file)

# Generar los histogramas y guardar los archivos temporales
for i in tqdm(range(num_steps)):
    save_histogram(a[0], a[1], i, t)

# Crear una animación a partir de las imágenes temporales
def animate(i):
    img = plt.imread(image_files[i])
    plt.imshow(img)
    plt.title(f'Time Step {i}/{num_steps}')
    plt.axis('off')

# Crear la animación
ani = FuncAnimation(plt.figure(), animate, frames=num_steps, interval=1000)
plt.show()
# Guardar la animación en un archivo de video
#ani.save('growingCoil.mp4', writer='ffmpeg', fps=200) 

# Limpiar los archivos temporales
cleanup_files(image_files)
