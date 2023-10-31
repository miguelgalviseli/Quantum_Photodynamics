from radiation import *
N=750
omega_c,omega_l = 0.5, 0.5
E0 = 0.05
ini = [ 0]  
num_steps = 6000
T=2000


average_RWA(N, omega_l,omega_c, E0, ini, num_steps,T)
average(N, omega_l,omega_c, E0, ini, num_steps,T)
