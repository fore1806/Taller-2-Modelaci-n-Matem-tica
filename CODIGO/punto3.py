# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 20:38:23 2022

@author: felip
"""

import numpy as np
import matplotlib.pyplot as plt


def ThomasSolve(a,b,c,d):
    n = len(d)
    
    # Modifica los coeficientes de la primera fila
    c[0] /= b[0]  # Posible división por cero
    d[0] /= b[0]
    
    #Creamos el espacio de respuesta
    x = []

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

def producto_matrices(a, b):
    filas_a = len(a)
    filas_b = len(b)
    columnas_a = len(a[0])
    columnas_b = len(b[0])
    if columnas_a != filas_b:
        return None
    # Asignar espacio al producto. Es decir, rellenar con "espacios vacíos"
    producto = []
    for i in range(filas_b):
        producto.append([])
        for j in range(columnas_b):
            producto[i].append(None)
    # Rellenar el producto
    for c in range(columnas_b):
        for i in range(filas_a):
            suma = 0
            for j in range(columnas_a):
                suma += a[i][j]*b[j][c]
            producto[i][c] = suma
    return producto
    
D = 1.8 #[m^2/s]
L = 0.9 #[m]
U = 0.2 #[m/s]
Cref = 0.1 #[kg/m^3]
k = 5 #[1/s]


n=10#Puntos espaciales de la malla
nT = 5000 #Puntos temporales

Ttotal = 10 #[s]


sizeA = (n,n)
sizeB = (n,n)
sizeC = (n,1)
sizeD = (n,1)

sizeY = (n,nT)
#sizeY = (nT,n)

DeltaX = L/n
DeltaT = Ttotal/nT

c1 = (1/DeltaT)+((U**2)/(2*DeltaX))+((L*U)/(2*D))+((L*k)/(2*U))+((D*U)/(L*(DeltaX**2)))
c2 = (D*U)/(L*(DeltaX**2))
c3 = (1/(2*DeltaX))+((D*U)/(2*L*(DeltaX**2)))
c4 = (1/DeltaT)+((k*L)/(2*U))+((D*U)/(L*(DeltaX**2)))
c5 =  (1/(2*DeltaX))-((D*U)/(2*L*(DeltaX**2)))
c6 =  (1/DeltaT)-((U**2)/(2*DeltaX))-((L*U)/(2*D))-((L*k)/(2*U))-((D*U)/(L*(DeltaX**2)))
c7 = (1/DeltaT)-((k*L)/(2*U))-((D*U)/(L*(DeltaX**2)))
c8 = ((U**2)/DeltaX)+(L*U/D)

A = np.zeros(sizeA) #Matriz A
B = np.zeros(sizeB) #Matriz B
C = np.zeros(sizeC) #Vector c
DV = np.zeros(sizeD) #Vector independiente

time = np.linspace(0, Ttotal, nT)
x = np.linspace(0, L, n)

a = np.zeros(sizeD)
b = np.zeros(sizeD)
c = np.zeros(sizeD)

Y = np.zeros(sizeY)
yN = np.zeros(sizeD) #Vector de concentración actual


#Llenado inicial matrices A y B y vector C
for i in range (n):
    arrayA=np.zeros(n)
    arrayB=np.zeros(n)
    if (i==0):
        arrayA[0] = c1
        arrayA[1] = -c2
        arrayB[0] = c6
        arrayB[1] = c2
        C[i] = c8
    elif (i==(n-1)):
        arrayA[n-1] = c4
        arrayA[n-2] = -c2
        arrayB[n-1] = c7
        arrayB[n-2] = c2
    else:
        arrayA[i-1]= -c3
        arrayA[i]= c4
        arrayA[i+1]= c5
        arrayB[i-1]= c3
        arrayB[i]= c7
        arrayB[i+1]= -c5
        C[i] = 0
    A[i] = arrayA
    B[i] = arrayB
    
    #Llenamos los vectores para el algoritmo de Thomas
    b[i] = arrayA[i]
    if (i>0):
        a[i] = arrayA[i-1]
    else:
        a[i] = 0
        
    if(i<(n-1)):
        c[i] = arrayA[i+1]
    else:
        C[i] = 0

for k in range (len(time)):
        productoo = producto_matrices(B,yN)
        DV = np.array(productoo) + C
        yOut = ThomasSolve(a,b,c,DV)
        for i in range (len(yOut)):
            if((k+1)<len(time)):
                Y[i,k+1] = yOut[i]
                yN[i] = yOut[i]


y10 = Y[:,nT-1]
y10 = np.multiply(Cref,y10)
plt.plot(x,y10)