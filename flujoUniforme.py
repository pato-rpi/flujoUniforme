import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# Inicializamos variables

Q = 1.20                # Caudal (m³/s)
n = 0.0120              # Coeficiente Manning 
b = 0.80                # Ancho basal (m)
z = 0.0                 # Talud Transversal (z:1=h:v) (m/m)
i = 0.20000             # Pendiente longitudinal (m/m)
l_inf = 0.1             # límite inferior de búsqueda de raíz.
l_sup = 10              # límite superior de búsqueda de raíz.

#Defino ecuación de manning trapezoidal.
def ecuacionAlturaNormal(h):
    return Q*n/(i**(0.5))-area(h)*radioHidraulico(h)**(2/3)

#Defino ecuación de altura critica trapezoidal.
def ecuacionAlturaCritica(h):
    return 1-Q**2*largoSuperficial(h)/(9.8*area(h)**3)

#Defino parámetros propios del escurrimiento
def area(h):
    return b*h + z*h*h

def perimetroMojado(h):
    return b + 2*h*((1)**2+(z)**2)**(0.5)

def radioHidraulico(h):
    return area(h)/perimetroMojado(h)

def largoSuperficial(h):
    return b + 2*z*h

def velocidad(h):
    return Q/area(h)

def froude(h):
    return (Q**2*largoSuperficial(h)/(9.8*area(h)**3))**(0.5)

def energiaEspecifica(h):
    return velocidad(h)**2/(2*9.8)+h

def momenta(h):
    return Q**2/(9.8*area(h))+centroGravedad(h)*area(h)

def centroGravedad(h):
    return 0
    
#Resolvemos la ecuación de manning obteniendo la raiz en función de un [a,b]
alturaNormal, info = brentq(ecuacionAlturaNormal, l_inf, l_sup, full_output=True)
print("Convergencia: ", info.converged)
print("Altura normal: ", alturaNormal, "m")
print("Área de escurrimiento: ", area(alturaNormal), "m²")
print("Perímetro mojado: ", perimetroMojado(alturaNormal), "m")
print("Radio hidráulico: ", radioHidraulico(alturaNormal), "m")
print("Ancho superficial: ", largoSuperficial(alturaNormal), "m")
print("Velocidad: ",velocidad(alturaNormal), "m/s")
print("Froude: ", froude(alturaNormal), "-")
print("Energía específica: ", energiaEspecifica(alturaNormal) , "m")
print("Momenta: ", "pendiente determinación Centro de Gravedad")
print("-----------------------------------------")

#Resolvemos la ecuación de manning obteniendo la raiz en función de un [a,b]
alturaCritica, info = brentq(ecuacionAlturaCritica, l_inf, l_sup, full_output=True)
print("Convergencia: ", info.converged)
print("Altura crítica: ", alturaCritica, "m")
print("Área de escurrimiento: ", area(alturaCritica), "m²")
print("Perímetro mojado: ", perimetroMojado(alturaCritica), "m")
print("Radio hidráulico: ", radioHidraulico(alturaCritica), "m")
print("Ancho superficial: ", largoSuperficial(alturaCritica), "m")
print("Velocidad: ",velocidad(alturaCritica), "m/s")
print("Froude: ", froude(alturaCritica), "-")
print("Energía específica: ", energiaEspecifica(alturaCritica) , "m")
print("Momenta: ", "pendiente determinación Centro de Gravedad")
