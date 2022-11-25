# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 23:10:10 2022

@author: felip
"""

import numpy as np
import matplotlib.pyplot as plt


def ThomasSolve(a,b,c,d):
    n = len(D)
    
    # Modifica los coeficientes de la primera fila
    c[0] /= b[0]  # Posible división por cero
    d[0] /= b[0]

    for i in range(1, n):
        ptemp = b[i] - (a[i] * c[i-1])
        c[i] /= ptemp
        d[i] = (d[i] - a[i] * d[i-1])/ptemp

       # Sustitución hacia atrás
    x = [0 for i in range(n)]
    x[-1] = d[-1]

    for i in range(-2, -n-1, -1):
        x[i] = d[i] - c[i] * x[i+1]
    
    x = np.array(x)
    
    return x
    

u1 = 0
un = 0
n=5

sizeA = (n,n)
sizeD = (n,1)

DeltaY = 2/n
c1 = 1/(DeltaY**2)
c2 = -2/(DeltaY**2)
c3 = -3

A = np.zeros(sizeA)
D = np.zeros(sizeD)
a = np.zeros(sizeD)
b = np.zeros(sizeD)
c = np.zeros(sizeD)


#Llenado inicial
for i in range (n):
    array=np.zeros(n)
    if (i==0):
        array[0] = 1
        D[0] = u1
        a[0] = 0;
        c[0] = 0;
    elif (i==(n-1)):
        array[n-1] = 1
        D[n-1] = un
        a[n-1] = 0
        c[n-1] = 0
    else:
        array[i-1]= c1
        array[i+1]= c1
        array[i]= c2
        D[i] = c3
        a[i] = c1
        c[i] = c1
    A[i] = array
    b[i] = array[i]


u = ThomasSolve(a,b,c,D)

x = 1
y = np.linspace(0,2,n)
xmesh,ymesh = np.meshgrid(x,y)

#flechitas
plt.figure(1)
plt.quiver(xmesh,ymesh,u,0)
# Curva
plt.figure(2)
plt.plot(u,y)
plt.show()

print(u)

#Procedimiento analitico: 
Re=300
k=-3/Re
miu=1 #viscosidad
um=2 #m/s
h=5 #m
uan=(y**2)*(k/(2*miu))+y*(0.25*um/(2*h)-(h*k/miu))
plt.figure(3)
plt.plot(uan,y)
plt.show
    
