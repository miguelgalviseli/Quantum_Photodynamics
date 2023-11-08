from population_average_const_env_retraso_negative import *
N=50
omega_c, omega_0, omega_l = 0.05, 0.05, 0.05
g = 0.01
E0 = 0.02
n=6
area = "inversion"
ini = ["e", 0]  
num_steps = 2000
tg= (n)*np.pi/(2*g)-200
tf=2*tg
ti=tg

#population_inversion_variedades(N, omega_l, omega_0, omega_c, g, E0, n, variedades, area, ini, num_steps)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
print(tf)
print(ti)
print(tg)
#population_inversion(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,tf,ti,tg)
population_inversion_all(N, omega_l, omega_0, omega_c,g, E0, n, area, ini, num_steps,tf,ti,tg)
#population_inversion_RWA_coherent(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps,retraso)
#average_photons_variedad(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)
#population_inversion_all(N, omega_l, omega_0, omega_c, g, E0, n, area, ini, num_steps)


