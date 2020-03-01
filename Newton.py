import numpy as np
import sympy as sp

funciones = [] #Guardar las funciones
variables = [] #Guardar las variables de las funciones
punto = [] #Guardar el punto inicial y posteriormente la solucion

n = int(input("Cuantas variables tiene tu sistema: "))

for i in range(n):
    funciones.append(input("función {}: ".format(i+1))) #Agregar n funciones a la lista

for i in range(n):                                           #
    for j in funciones[i]:                                   #
        if j <= 'z' and j >= 'a':                            #Agregar las variables a la lista
            variables.append(j)                              #y si alguna variable ya está, esta
            if j in variables and len(variables) > n:        #se elimina
                variables.pop(-1)                            #
variables.sort() #Ordenar las variables en orden alfabético

for i in range(n):
    punto.append(float(input("Ingresa el valor de {} en x0: ".format(variables[i])))) #Ingresa
                                                                                #los valores del
                                                                                #punto inicial

for i in variables:       #Hacer las
    i = sp.symbols(i)     #variables símbolos
