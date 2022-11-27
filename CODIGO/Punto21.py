# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 02:35:14 2022

@author: felip
"""

import numpy as np
import matplotlib.pyplot as plt

# Parámetros 
D=1.8  #m²/s    
L=0.9  #m      
U=0.2  #m/s    
Cref=0.1 #kg/m³ 
k=5 #s-1    


#paso en x
dx = 0.02 #m

#paso en t
tTotal = 10.0 #s
dt = 0.01 #s

#Puntos espaciales de la malla
n = int((L-0)/dx)+1

#Puntos temporales
nT = int((tTotal-0)/dt)+1

# vector de posiciones
xn = np.zeros(n)
xn[0] = 0
for i in range(1,n):
    xn[i] = xn[i-1]+dx


#Coeficientes de la matriz
C1=D/(U*L)
C2=(k*L)/U
C3=(C1*dt)/((dx**2))
C4=dt/(2*dx)
C5=C2*dt
C6=(C1**-1)*2*dx
C7=(1/dt)-(k*L/2*U)
C8=((U**2)/dx)+(L*U/D)


# MATRIZ A
A = np.zeros([n,n])

#primer nodo
A[0,0] = 1+(2*C3)+C5+(C6*(C3+C4))
A[0,1] = -2*C3

#Matriz Trigonal
for i in range(1,n-1):
    A[i,i-1] = -C3-C4
    A[i,i] = 1+(2*C3)+C5
    A[i,i+1] = -C3+C4
    
#ultimo nodo
A[-1,-2] = -2*C3
A[-1,-1] = 1+2*C3+C5

# MATRIZ Y
Y = np.zeros([n,nT])

# MATRIZ B
B = np.zeros([n,1])

# MATRIZ C
C = np.zeros([n,1])
C[0] = C6*(C3+C4)


# Solucionar Sitema AY=B+C
for i in range(1,nT):
    Y[:,[i]] = np.linalg.solve(A,B)+ np.linalg.solve(A,C)
    B = Y[:,[i]]

#dimensionalizar y* para obtener y.
Y=Y*Cref

#grafica en t=10s
plt.figure(1)
plt.plot(xn, Y[:,[14]],'b')
plt.title("Concentración de reactivo en instante t=10s")
plt.xlabel('Posición en reactor, x (m)')
plt.ylabel('Concentración, y (kg/m^3)')
plt.grid()


#grafica de ys en diferentes instantes de tiempo
plt.figure(2)
plt.plot(xn, Y[:,[14]],'b',label='t = 10s')
plt.plot(xn, Y[:,[10]],'r',label='t = 7.1s')
plt.plot(xn, Y[:,[6]],'g',label='t = 4.2s')
plt.plot(xn, Y[:,[3]],'black',label='t = 2.1s')
plt.plot(xn, Y[:,[2]],'orange',label='t = 1.4s')
plt.plot(xn, Y[:,[1]],'magenta',label='t = 0.7s')
plt.plot(xn, Y[:,[0]],'brown',label='t = 0s')
plt.legend(loc='upper right',ncol=2)
plt.title("Concentración de reactivo en diferentes instantes de tiempo")
plt.xlabel('Posición en reactor, x (m)')
plt.ylabel('Concentración, y (kg/m^3)')
plt.grid()
