import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# Inicializamos variables

Q = 30                  # Caudal (m³/s)
n = 0.0225              # Coeficiente Manning 
b = 8.20                # Ancho basal (m)
z = 0.5                 # Talud Transversal (z:1=h:v) (m/m)
i = 0.0005              # Pendiente longitudinal (m/m)
l_inf = 0               # límite inferior de búsqueda de raíz.
l_sup = 10              # límite superior de búsqueda de raíz.

#Defino ecuación de manning trapezoidal.
def f(h):
    A = b*h + z*h*h
    P = b + 2*h*((1)**2+(z)**2)**(0.5)
    R = A/P
    return Q*n/(i**(0.5))-A*R**(2/3)

#Resolvemos la ecuación de manning obteniendo la raiz en función de un [a,b]
zero, info = brentq(f, l_inf, l_sup, full_output=True)
print(zero, info.converged)

#Grafica la raiz encontrada.
h = np.linspace(l_inf, l_sup, 100)
plt.plot(h, f(h))
plt.axhline(color = 'k')
plt.axvline(x=zero, color = 'r', **{'linestyle': 'dashed'})
plt.show()
