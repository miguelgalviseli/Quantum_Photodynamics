from population_average_const_env import *
import matplotlib.pyplot as plt
N=35
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05

g = 0.001
E0 = 0.002
n=1
area = "inversion"
ini = ["e", 0]  
num_steps = 7000

variedades = [1,2,3,4,5,6,7,8,9,10,11,12]
tiempos = [5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900]
T = (2*n+1)*np.pi/(g)
#T=1000
#print(T)
#T=100
t = np.linspace(0,T,num_steps)
print(T)

#print(T)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tiempos)
#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)


