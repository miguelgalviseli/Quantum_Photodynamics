from population_average_sinsquare_env import *
N=50
omega_c, omega_0, omega_l = 0.5, 0.5, 0.5
g = 0.001
E0 = 0.002
n=1
area = "inversion"
ini = ["e", 0]  
num_steps = 12000
variedades = [1,2,3]

population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#print(len(a))