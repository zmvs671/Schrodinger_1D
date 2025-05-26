import numpy as np

def segundaDerivadaCentrada(f, x, h=1e-5):
    """
    f : Función a derivar, debe aceptar un float y devolver un float.
    x : Punto en el cual se evalúa la segunda derivada.
    h : Tamaño del paso para las diferencias finitas.
    """
    return (f(x + h) - 2 * f(x) + f(x - h)) / h**2

def numerov(V, x_min, x_max, h, E, psi_0=0.0, psi_1=1e-5):
    """
    - V: función del potencial V(x)
    - x_min, x_max: intervalo de posiciones
    - h: paso de integración
    - E: energía tentativa
    - psi_0, psi_1: condiciones iniciales (\psi(x0), \psi(x₀ + h))
    """
    m = 1.0     # masa reducida 
    hbar = 1.0  # Constante de Planck reducida 

    # Número de pasos espaciales a partir del intervalo y el tamaño del paso h
    N = int((x_max - x_min) / h)

    # Genera un array de N puntos espaciados uniformemente entre x_min y x_max
    x = np.linspace(x_min, x_max, N)

    # Calcula k^2(x) = (2m/\hbar^2)(E - V(x)) en cada punto del dominio
    k2 = 2 * m * (E - V(x)) / hbar**2

    psi = np.zeros(N)

    # Condiciones iniciales: \psi(x0) = psi_0 y \psi(x0 + h) = psi_1
    psi[0] = psi_0
    psi[1] = psi_1

    # Iteración del método de Numerov para obtener \psi(x) en todo el dominio
    for n in range(1, N - 1):
        psi[n + 1] = (
            (2 * (1 - (5 * h**2 * k2[n] / 12)) * psi[n] -  # Término central
            (1 + (h**2 * k2[n - 1] / 12)) * psi[n - 1]) / # Término anterior
            (1 + (h**2 * k2[n + 1] / 12))                  # Término siguiente (denominador)
        )
    return x, psi #x es un array de posiciones, psi es un array de eigenfunciones

def rk4d(f, x0, a, b, h=0.01):
  '''
  f: función del sistema diferencial dx/dt = f(x, t)
  x0: condición inicial (valor de x en t = a)
  a, b: intervalo en t 
  h: paso de integración 
  '''
  xlista = []  # Lista para guardar posiciones
  tlista = np.arange(a, b + h, h)  # Lista de puntos en el tiempo desde a hasta b con paso h

  x = x0  # Establecemos la condicion inicial

  for i in tlista:
    xlista.append(x) 

    # Método de Runge-Kutta 4to orden:
    k1 = h * f(x, i)                            # Primera estimacion (como Euler)
    k2 = h * f(x + (k1 / 2), i + h / 2)         
    k3 = h * f(x + (k2 / 2), i + h / 2)         
    k4 = h * f(x + k3, i + h)                   # Evaluación al final del paso con k3

    # Combinacion ponderada para actualizar posicion
    x = x + (k1 + 2*k2 + 2*k3 + k4) / 6

  return tlista, xlista  # Regresamos la lista de tiempos y la solucion x(t)

def secante(f, x0, x1, epsilon=1e-6):
  '''
  f: función de la que queremos encontrar una raíz
  x0, x1: primeras dos aproximaciones iniciales (se requieren dos)
  epsilon: tolerancia para determinar cuándo f(x) es suficientemente cercana a 0
  '''
  f0 = f(x0)  # Evaluamos la funcion en el primer punto
  f1 = f(x1)  # Evaluamos la funcion en el segundo punto

  while True:
    # Para evitar divisiones inseguras o errores de precision, usamos **(-1) en vez de "/"
    x2 = x1 - f1 * (x1 - x0) * (f1 - f0) ** (-1)

    f2 = f(x2)  # Evaluamos la funcion en la nueva aproximacion

    if abs(f2) < epsilon:
      break  # Criterio de convergencia
    
    # Actualizamos las variables para la siguiente iteracion
    f0 = f1
    f1 = f2
    x0 = x1
    x1 = x2

    return x2  # Regresa la mejor aproximacion de la raiz

def QRdescomp(A):
    '''
    Descompone una matriz cuadrada A en el producto QR usando Gram-Schmidt.
    A : ndarray (n x n) Matriz real cuadrada a descomponer.
    '''
    n = A.shape[0]  # Se obtiene el número de filas = columnas de A
    
    # Inicializamos Q y R como matrices de ceros
    Q = np.zeros_like(A, dtype=float)
    R = np.zeros((n, n), dtype=float)

    # Gram-Schmidt modificado
    for i in range(n):
        u = A[:, i].copy()  # Tomamos la i-ésima columna de A

        # Restamos las proyecciones sobre los vectores ortonormales anteriores
        for j in range(i):
            R[j, i] = np.dot(Q[:, j], A[:, i])  # Producto interno q_j · a_i
            u -= R[j, i] * Q[:, j]              # Restamos la componente en dirección de q_j

        R[i, i] = np.linalg.norm(u)  # Norma del vector ortogonal resultante
        Q[:, i] = u / R[i, i]        # Normalizamos para obtener q_i

    return Q, R  # Se devuelve la matriz ortogonal Q y la triangular superior R

def QR(A, tol=1e-6):
    '''
    Calcula los eigenvalores y eigenvectores de una matriz simétrica real usando el algoritmo QR.
    A : ndarray (n x n). Matriz real y simétrica cuyas raíces (eigenvalores) se desean calcular.
    tol : Tolerancia para determinar la convergencia. El algoritmo termina cuando 
        todos los elementos fuera de la diagonal son menores que esta tolerancia.
    '''
    # Obtenemos el tamaño de la matriz cuadrada A
    n = A.shape[0]
    
    # Inicializamos Qtotal como la matriz identidad de tamaño n x n
    Qtotal = np.eye(n)
    
    while True:
        # Calculamos la descomposición QR de la matriz Ak
        Q, R = QRdescomp(A)
        
        # Obtenemos la nueva matriz Ak = R @ Q (producto de R por Q)
        Ak = R @ Q
        
        # Acumular el producto de Qs para construir la matriz de eigenvectores
        Qtotal = Qtotal @ Q
        
        # Construimos la matriz que contiene solo los elementos fuera de la diagonal
        offDiagonal = Ak - np.diag(np.diagonal(Ak))
        
        # Criterio de convergencia: los elementos fuera de la diagonal son suficientemente pequeños
        if np.all(np.abs(offDiagonal) < tol):
            break  

    # Los eigenvalores están en la diagonal de la matriz Ak final
    eigenvalores = np.diagonal(Ak)
    
    # Los eigenvectores se guardan en las columnas de Qtotal
    eigenvectores = Qtotal
    
    return eigenvalores, eigenvectores