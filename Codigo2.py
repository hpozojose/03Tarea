'''
Este script contiene la resolucion de la parte 2 de la tarea 3, donde se pide
integrar la ecuacion de lorentz mediante runge kutta
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode

def f_to_solve(t, V):
    return [si*(V[1]-V[0]), V[0]*(rho-V[2])-V[1],V[0]*V[1]-beta*V[2]]


# Condiciones iniciales
t0 = 1e-3
V0 = [1.22, 1.22, 1.22]
si=10
rho=28
beta=8/3

# creamos el 'resolvedor'
r = ode(f_to_solve)
r.set_integrator('dopri5')
# r.set_integrator('dopri5')
r.set_initial_value(V0)

# guardamos las variables a medida que progresamos
t = 10000
t_values = np.linspace(t0, 40 * np.pi, t)
x_values = np.zeros(t)
y_values = np.zeros(t)
z_values = np.zeros(t)

for i in range(len(t_values)):
    r.integrate(t_values[i])
    x_values[i], y_values[i], z_values[i] = r.y

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

ax.plot(x_values, y_values, z_values)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
