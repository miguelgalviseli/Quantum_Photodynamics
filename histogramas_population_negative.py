
import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import cv2
import os
from population_average_const_env_retraso_negative import *
N=50
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=10
area = "inversion"
ini = ["e", 0]  
num_steps = 1000
tg= (2*n+1)*np.pi/(2*g)-2000
tf=2*tg
ti=tg

#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tg)
print(ti)
print(tf)

#print(T)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
t = np.linspace(0,max([tf,tg+ti]),num_steps)

a = population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg)
print(len(a[2]))

#La idea es graficar los histogramas de los datos de a para cada tiempo t y luego hacer una animación
#Defino la función que va a hacer el histograma para cada tiempo
def histograma_tiempo(e, g,intensidad, tiempo,t):
    #Debo tomar el enesimo elemento de cada sublista de la lista de listas
    #Para eso, voy a hacer un ciclo for que recorra cada sublista, y que tome el enesimo elemento de cada sublista


    #Defino la lista donde voy a guardar los enesimos elementos de cada sublista
    lista_enelementose = []
    lista_enelementosg = []

    #Defino el ciclo for que recorre cada sublista
    for i in range(len(e)):
        #Defino el enesimo elemento de cada sublista
        enelementoe = e[i][tiempo]
        enelementog = g[i][tiempo]

        #Agrego el enesimo elemento a la lista de enelementos
        lista_enelementose.append(enelementoe)
        lista_enelementosg.append(enelementog)

    #Defino la lista de n
    lista_n = []

    #Defino el ciclo for que va a definir la lista de n
    for i in range(len(e)):
        #Defino n
        n = i
        #Agrego n a la lista de n
        lista_n.append(n)

    #Defino el histograma
    
    
    plt.figure(figsize=(10,6)) # Crear una nueva figura en cada iteración
    #Pintemos un 1 en el eje para todos los histogramas
    
    plt.bar(lista_n, lista_enelementose, color = '#173F5F', edgecolor = 'black', width = 1)
    plt.bar(lista_n, lista_enelementosg, color = 'orange', edgecolor = 'black', width = 1)
    plt.xlabel('Variedad n',fontsize=15)
    plt.ylabel('Inversión',fontsize=15)
    plt.title('Jaynes-Cummings y radiación. t = ' + str(np.round(t,2)),fontsize=15)
    #plt.text(N-5, 0.9, "Intensidad = {}".format(round(intensidad[tiempo],6)), fontsize=10, ha='center')
    #plt.annotate("Texto en la parte superior", xy=(t, intensidad[tiempo]), xycoords='axes fraction', fontsize=12, ha='center')
    plt.bar(N+1, intensidad[tiempo]*1000, color = 'red', edgecolor = 'black', width = 1, label='Intensidad')
    plt.axhline(y=1, color='g', linestyle='-')
    plt.axhline(y=-1, color='g', linestyle='-')
    plt.legend(fontsize=6)
    #guardar la imagen en una carpeta 
    
    #Quiero que queden guardadas en una carpeta distinta a la del código
    plt.savefig('Histogramas_f1/' + 'histograma_tiempo_' + str(tiempo) + '.png')
    plt.close()

#Defino el ciclo for que va a hacer el histograma para cada tiempo
for i in tqdm(range(len(a[0][0]))):
    histograma_tiempo(a[0],a[1],a[2], i, t[i])

# Ruta donde se encuentran las imágenes
images_folder = 'Histogramas_f1/' 

# Nombre del video de salida
output_video = "video_salida2.mp4"

# Obtener la lista de nombres de archivos en orden

image_names = sorted(
    [img for img in os.listdir(images_folder) if img.endswith(".png")],
    key=lambda x: int(x.split("_")[2].split(".")[0])
)


# Leer la primera imagen para obtener sus dimensiones
first_image = cv2.imread(os.path.join(images_folder, image_names[0]))
height, width, layers = first_image.shape

# Configurar el codec y crear un objeto VideoWriter
fourcc = cv2.VideoWriter_fourcc(*"mp4v")  # Puedes cambiar el codec según tus necesidades
#Quiero que el video duré menos tiempo
video = cv2.VideoWriter(output_video, fourcc, 20, (width, height))

# Agregar cada imagen al video
for image_name in image_names:
    img_path = os.path.join(images_folder, image_name)
    img = cv2.imread(img_path)
    video.write(img)


# Liberar recursos
video.release()
cv2.destroyAllWindows()

print("Video creado exitosamente:", output_video)