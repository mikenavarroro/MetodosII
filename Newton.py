import numpy as np
import sympy as sp

funciones = [] #Guardar las funciones
variables = [] #Guardar las variables de las funciones

#Hacer símbolos el abecedario para poder derivar
a, b, d, f, j, k = sp.symbols('a b d f j k')
m, q, r, u, v, w, x, y, z = sp.symbols('m q r u v w x y z')

abecedario = [a, b, d, f, j, k, m, q, r, u, v, w, x, y, z]

#Pedir el número de variables del sistema
n_var = int(input("Cuantas variables tiene tu sistema: "))

punto = [] #Guardar el punto inicial y posteriormente la solución

for i_var in range(n_var):
    funciones.append(input("función {}: ".format(i_var+1))) #Agregar n funciones a la lista

for i_var in range(n_var):
    for j_var in funciones[i_var]:                                   #
        if j_var >= 'A' and j_var <= 'Z':                            #Convertir las variables
            j_var = j_var.lower()                                    #en minúsculas
        if j_var >= 'a' and j_var <= 'z':                            #Agregar las variables a la lista
            variables.append(j_var)                                  #y si alguna variable ya está, esta
            if j_var in variables and len(variables) > n_var:        #se elimina
                variables.pop(-1)                                    #
variables.sort() #Ordenar las variables en orden alfabético

#Ingresar el punto inicial
for i_var in range(n_var):
    punto.append(float(input("Ingresa el valor de {} en x0: ".format(variables[i_var]))))

punto = np.asarray(punto)

#Pedir la tolerancia
tol = float(input("Tolerancia: "))
#Número máximo de iteraciones
iter = int(input("Numero máximo de iteraciones: "))

Jacobi = np.empty((n_var, n_var)) #Matriz jacobina vacía
func = [] #Lista para poder derivar
der = [] #Lista que guarda las derivadas
no_iteraciones = 1 #Número de iteraciones hechas
error = 1 #Error
puntos_anteriores = [] #Cálculo del error

#Hacer símbolos las funciones
for i_var in range(n_var):
    func.append(sp.sympify(funciones[i_var]))

#Derivar las funciones con respecto a las variables
for i_var in range(n_var):
    for j_var in range(n_var):
        der.append(func[i_var].diff(variables[j_var]))

var_abc = [] #Guardar las variables del abecedario que se ocupan en símbolos

#Agregar las variables usadas a la lista var_abc para poder evaluar las funciones
for i_var in range(n_var):
    for j_var in variables[i_var]:
        var_abc.append(j_var)
        if j_var in variables and len(var_abc) > n_var:
            var_abc.pop(-1)
var_abc.sort() #Ordenar alfabéticamente las variables que se van a usar

#Inicio del método de Newton
while no_iteraciones <= iter:# or error > tol:
    #Evaluar la matriz jacobiana en el punto
    m_var = 0
    for i_var in range(n_var):
        for j_var in range(n_var):
            Jacobi[i_var, j_var] = sp.sympify(der[m_var]).subs(var_abc[j_var], punto[j_var])
            m_var += 1

    JacobiInv = np.linalg.inv(Jacobi) #Inversa de la matriz
    F = [] #Crear la matriz de F(Xk)

    #Evaluar las funciones en el punto
    for i_var in range(n_var):
        for j_var in range(n_var):
            if var_abc[j_var] in funciones[i_var]:
                F.append(sp.sympify(funciones[i_var]).subs(var_abc[j_var], punto[j_var]))

    Y = JacobiInv @ F

    puntos_anteriores.append(punto) #Guardar el punto en otra lista para el cálculo del error
    punto = punto - Y #Calcular X(k+1)

    #Calcular el error
    #if no_iteraciones != 1:
    #    p_var = 0
    #    q_var = 0
    #    while p_var <= len(punto):
    #        error = max(abs(punto[p_var] - puntos_anteriores[p_var, q_var]))
    #        p_var += 1

    print(no_iteraciones)
    no_iteraciones += 1

#Imprimir el resultado
if no_iteraciones > iter:
    print("El número de iteraciones fue excedido.")
    print("La solución es ", punto)
elif error < tol:
    print("Se ha alcanzado el límite de tolerancia establecida.")
    print("La solución es ", punto)
