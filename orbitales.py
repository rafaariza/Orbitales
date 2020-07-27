# Importamos todos los módulos necesarios

from numpy import sin, cos, pi, sqrt, shape, linspace, meshgrid, zeros, array
from math import factorial
import matplotlib.pyplot as plt

# Definimos funciones necesarias para nuestro cometido

def esfericas2cartesianas(r,theta,phi):

    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    
    return x, y, z

def legendre_polinomio(l,m,x):
    pmm = 1.0
    if m > 0:
        sign = 1.0 if m % 2 == 0 else -1.0
        pmm = sign*pow(factorial(2*m-1)*(1.0-x*x),((m/2)))

    if l == m:
        return pmm

    pmm1 = x*(2*m+1)*pmm
    if l == m+1:
        return pmm1

    for n in range(m+2,l+1):
        pmn = (x*(2*n-1)*pmm1-(n+m-1)*pmm)/(n-m)
        pmm = pmm1
        pmm1 = pmn

    return pmm1

# Ya tenemos todas las funciones definidas. Let's code!
    
l = int(input('Número cuántico azimutal l: '))

if l < 0:
    print('l no puede ser negativo')

m = int(input('Número cuántico magnético m: '))

if m>l or m<(-l):
    print('Número cuántico m pertenece a [-l,l]')

A = sqrt(((2*l+1)*factorial(l-abs(m)))/(4*pi*factorial(l+abs(m)))) # Constante de normalización

phi = linspace(0,2*pi,181)
theta = linspace(0,pi,91)

Phi,Theta = array(meshgrid(phi,theta))

if m>0:
        
    rho = pow(abs(sqrt(2)*A*cos(m*Phi)*legendre_polinomio(l,m,cos(Theta))),2)
        
elif m < 0:
        
    rho = pow(abs(sqrt(2)*A*sin(abs(m)*Phi)*legendre_polinomio(l,abs(m),cos(Theta))),2)
        
else:
        
    rho = pow(abs(A*legendre_polinomio(l,0,cos(Theta))),2)

x,y,z = esfericas2cartesianas(abs(rho),Theta,Phi)

# Figuras YAY!!

fig = plt.figure(f'Armónico {l},{m}')
ax = fig.add_subplot(111,projection='3d')

ax.plot_surface(x,y,z, cmap = 'Spectral')

plt.title(f'Armónico {l},{m}', fontsize = 14, fontweight = 'bold')
ax.set_xlabel('X', fontsize = 12, fontweight = 'bold')
ax.set_ylabel('Y', fontsize = 12, fontweight = 'bold')
ax.set_zlabel('Z', fontsize = 12, fontweight = 'bold')

plt.show()
