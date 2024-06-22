#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>
#include <gmpxx.h>
#include <vector>

using namespace std;

struct Punto {
    mpz_class x;
    mpz_class y;
};

/**
 * @file Utils.h
 * @brief Contiene diversas funciones de utilidad para operaciones matemáticas y criptográficas.
 */

/**
 * Genera un número aleatorio corto.
 * @param state Estado del generador de números aleatorios GMP.
 * @return Un número aleatorio de 64 bits.
 */
unsigned long int generate_random_limb(gmp_randstate_t state);

/**
 * Comprueba si un número es pseudoprimo fuerte en cierta base.
 * @param p Número entero a analizar.
 * @param base Base donde estudiar la pseudoprimalidad de p.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return Verdadero si p es pseudoprimo fuerte en base, falso en otro caso.
 */
bool strongPseudoprime(const mpz_class p, const mpz_class base, bool debug = false);

/**
 * Aplica el test de Miller Rabin sobre un primo cierto número de veces.
 * @param p Número entero a analizar.
 * @param it Número de iteraciones que aplicar la prueba.
 * @param state Estado del generador de números aleatorios GMP.
 * @return Verdadero si p es un primo con la confianza indicada, falso en otro caso.
 */
bool millerRabin(const mpz_class p, int it, gmp_randstate_t state);

/**
 * Función para generar un número primo de 'bits' bits.
 * @param bits Número de bits deseado del primo a generar.
 * @param state Estado del generador de números aleatorios GMP.
 * @return Un primo del tamaño especificado.
 */
mpz_class generate_prime(unsigned int bits, gmp_randstate_t state);

/**
 * Función para generar un número primo robusto de 'bits' bits.
 * @param bits Número de bits deseado del primo a generar.
 * @param state Estado del generador de números aleatorios GMP.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return Número de tipo mpz_class que representa el primo generado.
 */
mpz_class generate_strong_prime(unsigned int bits, gmp_randstate_t state, bool debug = false);

/**
 * Realiza la suma de dos puntos en una curva elíptica sobre Z_p.
 * @param s1 Primer punto a sumar.
 * @param s2 Segundo punto a sumar.
 * @param a Coeficiente a de la curva elíptica.
 * @param b Coeficiente b de la curva elíptica.
 * @param p Módulo primo de la curva elíptica.
 * @param inv Inverso multiplicativo.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return Punto resultante de la suma.
 */
Punto sumaCurvaEliptica(const Punto s1, const Punto s2, const mpz_class a, 
                        const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);

/**
 * Realiza la multiplicación de un punto por un escalar en una curva elíptica sobre Z_p.
 * @param f1 Punto a multiplicar.
 * @param f2 Escalar por el cual multiplicar el punto.
 * @param a Coeficiente a de la curva elíptica.
 * @param b Coeficiente b de la curva elíptica.
 * @param p Módulo primo de la curva elíptica.
 * @param inv Inverso multiplicativo.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return Punto resultante de la multiplicación.
 */
Punto multiplicacionCurvaEliptica(const Punto f1, const mpz_class f2, const mpz_class a, 
                                  const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);

/**
 * Calcula la potencia máxima k de un número dado como potencia suave.
 * @param k Número entero para calcular su máxima potencia suave.
 * @return La potencia máxima k como potencia suave.
 */
mpz_class maxKPotenciaSuave(const mpz_class k);

/**
 * Calcula la raíz cuadrada modular de 'a' módulo 'n'.
 * @param a Número entero a encontrar la raíz cuadrada.
 * @param n Módulo para la operación.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return La raíz cuadrada modular de 'a' módulo 'n'.
 */
mpz_class sqrt_mod(const mpz_class &a, const mpz_class &n, bool debug = false);

/**
 * Transpone una matriz representada como un vector de vectores de booleanos.
 * @param matrix Matriz a transponer.
 * @return La matriz transpuesta.
 */
vector<vector<bool>> transpose(vector<vector<bool>>& matrix);

/**
 * Realiza la eliminación gaussiana sobre una matriz de booleanos.
 * @param matrix Matriz sobre la cual aplicar la eliminación.
 * @return La matriz resultante después de la eliminación gaussiana.
 */
vector<vector<bool>> gaussian_elimination(const vector<vector<bool>> &matrix);

/**
 * Encuentra soluciones para una matriz dada de ecuaciones booleanas.
 * @param matrix Matriz de ecuaciones booleanas.
 * @return Vectores de soluciones para la matriz dada.
 */
vector<vector<bool>> find_solutions(const vector<vector<bool>>& matrix);

/**
 * Calcula los cocientes de la fracción continua para la relación num/den.
 * @param num Numerador de la fracción.
 * @param den Denominador de la fracción.
 * @return Vector de mpz_class con los cocientes de la fracción continua.
 */
vector<mpz_class> cocientes_fraccion_continua(const mpz_class &num, const mpz_class &den);

/**
 * Resuelve una ecuación cuadrática de la forma ax^2 + bx + c = 0.
 * @param p Primera solución de la ecuación cuadrática.
 * @param q Segunda solución de la ecuación cuadrática.
 * @param a Coeficiente a de la ecuación cuadrática.
 * @param b Coeficiente b de la ecuación cuadrática.
 * @param c Término constante de la ecuación cuadrática.
 * @return Verdadero si la ecuación tiene soluciones, falso en otro caso.
 */
bool resuelveEcuacionCuadratica(mpz_class &p, mpz_class &q, mpz_class a, mpz_class b, mpz_class c);

#endif // UTILS_H


