import numpy as np
import sympy as sp

variables = [] #Guardar las variables de las funciones

#Hacer símbolos el abecedario para poder derivar
a, b, c, d, e, f, g, h, i, j, k = sp.symbols('a b c d e f g h i j k')
l, m, n, o, p, q, r, s, t, u, v = sp.symbols('l m n o p q r s t u v')
w, x, y, z = sp.symbols('w x y z')

abecedario = [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z]

#Pedir el número de variables del sistema
n_var = int(input("Cuantas variables tiene tu sistema: "))

funciones = np.array((n_var, 1)) #Guardar las funciones
punto = np.array((n_var, 1)) #Guardar el punto inicial y posteriormente la solución

for i_var in range(n_var):
    funciones[i_var, 0] = input("función {}: ".format(i_var+1)) #Agregar n funciones a la lista

for i_var in range(n_var):                                           #
    for j_var in funciones[i_var, 0]:                                #
        if j_var >= 'A' and j_var <= 'Z':                            #Convertir las variables
            j_var = j_var.lower()                                    #en minusculas
        if j_var >= 'a' and j_var <= 'z':                            #Agregar las variables a la lista
            variables.append(j_var)                                  #y si alguna variable ya está, esta
            if j_var in variables and len(variables) > n_var:        #se elimina
                variables.pop(-1)                                    #
variables.sort() #Ordenar las variables en orden alfabético

#Ingresar el punto inicial
for i_var in range(n_var):
    punto[i_var, 0] = float(input("Ingresa el valor de {} en x0: ".format(variables[i_var])))

#Pedir la tolerancia
tol = float(input("Tolerancia: "))
#Número máximo de iteraciones
iter = int(input("Numero máximo de iteraciones: "))

Jacobi = np.empty((n_var, n_var)) #Matriz jacobina vacía
func = [] #Lista para poder derivar
der = [] #Lista que guarda las derivadas
o_var = 1 #Número de iteraciones hechas
error = 1 #Error
puntos_anteriores = [] #Cálculo del error

#Hacer símbolos las funciones
for i_var in range(n_var):
    func.append(sp.sympify(funciones[i_var, 0]))

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
while o_var < iter or error > tol:
    #Evaluar la matriz jacobiana en el punto
    m_var = 0
    for i_var in range(n_var):
        for j_var in range(n_var):
            Jacobi[i_var, j_var] = sp.sympify(der[m_var]).subs(var_abc[j_var], punto[j_var, 0])
            m_var += 1

    JacobiInv = np.linalg.inv(Jacobi) #Inversa de la matriz

    F = np.empty((n_var, 1)) #Crear la matriz de F(Xk)

    #Evaluar las funciones en el punto
    l_var = 0
    for i_var in range(n_var):
        for j_var in range(n_var):
            F[i_var, 1] = sp.sympify(funciones[l_var, 0]).subs(var_abc[j_var], punto[j_var, 0])
            l_var += 1

    Y = JacobiInv @ F #Multiplicación de la inversa del jacobiano con el vector F(Xk)

    puntos_anteriores.append(punto) #Guardar el punto en otra lista para el cálculo del error

    punto = punto - Y #Calcular X(k+1)

    #Calcular el error
    p_var = 0
    for i_var in puntos_anteriores:
        error = max(abs(punto[p_var, 0] - i_var))

    o_var += 1

#Imprimir el resultado
if o_var > iter:
    print("El número de iteraciones fue excedido.")
    print("La solución es ", punto)
elif error < tol:
    print("Se ha alcanzado el límite de tolerancia establecida.")
    print("La solución es ", punto)
