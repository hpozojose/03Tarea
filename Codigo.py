
'''
Este script integra la una ecuacion de movimiento para el
oscilador de van der pol usando el metodo de Runge-Kutta
de tercer orden
'''

import numpy as np
import matplotlib.pyplot as plt
#Maria Jose Hernandez Pozo
#18.934.376-0


def f(y, dy):
    #funcion de van der pool simplificada y en matriz
    return dy, -y -mu*((y**2)-1)* dy

def get_k1(y_n, dy_n, h, f):
    #Definicion de k1
    f_eval = f(y_n, dy_n)
    return h * f_eval[0], h * f_eval[1]

def get_k2(y_n, dy_n, h, f):
    #definicon de k2
    k1 = get_k1(y_n, dy_n, h, f)
    f_eval = f(y_n + k1[0]/2, dy_n + k1[1]/2)
    return h * f_eval[0], h * f_eval[1]

def get_k3(y_n, dy_n, h, f):
    #Definicion de k3
    k1 = get_k1(y_n, dy_n, h, f)
    k2 = get_k2(y_n, dy_n, h, f)
    f_eval=   f(y_n-k1[0]- 2*k2[0], dy_n-k1[1]-2*k2[1])
    return h * f_eval[0], h * f_eval[1]

def rk3_step(y_n, dy_n, h, f):
    #Funcion que entrega el y'_(n+1/2)
    k1 = get_k1(y_n, dy_n, h, f)
    k2 = get_k2(y_n, dy_n, h, f)
    k3 = get_k3(y_n, dy_n, h, f)
    y_n1 = y_n + (k1[0] + 4 * k2[0] + k3[0])/6
    dy_n1 = dy_n + (k1[1] + 4 * k2[1] + k3[1])/6
    return y_n1, dy_n1

#definicion de salto o paso h, depende de el periodo de la funcion
N_steps = 40000
h = np.pi*20. / N_steps
y = np.zeros(N_steps)
dy = np.zeros(N_steps)
mu = 1.376

y[0] = 0.1
dy[0] = 0
#Guarda terminos en lista y y dy
for i in range(1, N_steps):
    y[i], dy[i] = rk3_step(y[i-1], dy[i-1], h, f)
#variable tiempo
t_rk = [h * i for i in range(N_steps)]

#graficos
plt.figure(1)
plt.clf()
plt.plot(t_rk, y, 'g')
plt.xlabel('variable s')
plt.ylabel('$y(s)$', fontsize=18)
plt.title('Grafico y(s) vs s')
plt.show()

plt.figure(2)
plt.clf()
plt.plot(y, dy, 'b')
plt.xlabel('variable y')
plt.ylabel('$dy(s)/ds$', fontsize=18)
plt.title('Grafico dy(s)/ds vs y')
plt.show()

y[0] = 4
dy[0] = 0
for i in range(1, N_steps):
    y[i], dy[i] = rk3_step(y[i-1], dy[i-1], h, f)


plt.figure(3)
plt.clf()
plt.plot(t_rk, y, 'b')
plt.xlabel('variable s')
plt.ylabel('$y(s)$', fontsize=18)
plt.title('Grafico y(s) vs s')
plt.show()

plt.figure(4)
plt.clf()
plt.plot(y, dy, 'b')
plt.xlabel('variable y')
plt.ylabel('$dy(s)/ds$', fontsize=18)
plt.title('Grafico dy(s)/ds vs y')
plt.show()
