# Soluciones Numéricas a la Ecuación de Schrödinger en 1D

Este proyecto implementa diversos métodos numéricos para resolver la ecuación de Schrödinger independiente del tiempo en una dimensión. A lo largo de este, se resuelven problemas en búsqueda de funciones de onda y niveles de energía en diferentes potenciales clásicos.

La ecuación de Schrödinger es fundamental en la mecánica cuántica, pero solo puede resolverse analíticamente en un número limitado de casos. Este repositorio plantea la construcción de algunos métodos que son puestos a prueba sobre problemas de los que se conoce su solución exacta y también se exploran sistemas donde las soluciones exactas no son tan prácticas.

Se aborda la ecuación estacionaria:

$$
\frac{d^2\psi(x)}{dx^2} + \frac{2m}{\hbar^2}(E - V(x))\psi(x) = 0
$$

## Métodos numéricos implementados

- **Diferencias finitas centradas**  
  Discretización directa de la ecuación para reformularla como un problema de autovalores.
  
- **Algoritmo QR**  
  Para la diagonalización de matrices tridiagonales en problemas de valores y vectores propios.

- **Método de Numerov**  
  Método iterativo de alta precisión para ecuaciones diferenciales del tipo \( y'' + k(x) y = 0 \).

- **Método de disparo (shooting)**  
  Ajuste de condiciones iniciales para encontrar autovalores mediante integración.

- **Método de Runge-Kutta de orden 4 (RK4)**  
  Método que permite resolver ecuaciones diferenciales ordinarias (EDO) de primer orden. La ecuación de Schrödinger si bien es de segundo grado, puede reformularse como un sistema de EDO de primer orden.


## Estructura del repositorio
├── metodos/ # Implementaciones de los métodos numéricos
├── problemas/ # Ejemplos aplicando los métodos a diferentes potenciales
├── requeriments.txt # Dependencias necesarias 
├── README.md # Este archivo

## Ejemplos incluídos
Todos los ejemplos están resueltos en notebooks individuales
- Potencial armónico --> armonico_Numerov.ipynb
- Verificación de eigenfunciones para el pozo infinito --> difinitas.ipynb
- Potencial cuadrado con paredes infinitas --> edoBase_Shooting.ipynb
- Pozo de ancho L con potencial no constante qr.ipynb

## Instala los requisitos 
```bash
pip install -r requeriments.txt
```
## Referencias

[1] Newman, M. E. J. (2013). *Computational physics*. CreateSpace Independent Publishing Platform.

[2] Burden, R. L., & Faires, J. D. (2010). *Numerical analysis* (9.ª ed.). Brooks/Cole, Cengage Learning.

[3] Caruso, F., Oguri, V., & Silveira, F. (2022). Applications of the Numerov method to simple quantum systems using Python. Revista Brasileira de Ensino de Física, 44, e20220098. https://doi.org/10.1590/1806-9126-RBEF-2022-0098

[4] Evans, M. (s.f.). Schroedinger Equation in Harmonic Potential [Código fuente]. Recuperado de https://mtdevans.com/js/sch.py.txt

## Créditos
Desarrollado por Silvia Maldonado como parte de un trabajo final para la asignatura de Mecánica Cuántica.

