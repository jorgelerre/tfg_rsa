#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>
#include <gmpxx.h>
#include <vector>

using namespace std;

struct Punto{
	mpz_class x;
	mpz_class y;
};

/*
 * Genera un numero aleatorio corto
 * @it Numero de iteraciones que aplicar la prueba
 * @return Un numero aleatorio de 64 bits
*/
unsigned long int generate_random_limb(gmp_randstate_t state);

/*
 * Comprueba si un numero es pseudoprimo fuerte en cierta base
 * @p Numero entero a analizar
 * @base Base donde estudiar la pseudoprimalidad de p
 * @return True si p es pseudoprimo fuerte en base, False en otro caso.
 */
bool strongPseudoprime(const mpz_class p, const mpz_class base, bool debug = false);

/*
 * Aplica el test de Miller Rabin sobre un primo cierto numero de veces
 * @p Numero entero a analizar
 * @it Numero de iteraciones que aplicar la prueba
 * @state El estado del generador de numeros aleatorios
 * @return True si p es un primo con confianza 1-4^{-it}, false en otro caso.
*/
bool millerRabin(const mpz_class p, int it, gmp_randstate_t state);

/**
 * Función para generar un número primo de 'bits' bits
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 * @return Un primo del tamanio especificado
 */
mpz_class generate_prime(unsigned int bits, gmp_randstate_t state);

/**
 * Función para generar un número primo robusto de 'bits' bits
 * 
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 * @return Dato de tipo mpz_class con el primo generado
 */
mpz_class generate_strong_prime(unsigned int bits, gmp_randstate_t state, bool debug = false);


//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
Punto sumaCurvaEliptica(const Punto s1, const Punto s2, const mpz_class a, 
					    const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);
					    
Punto multiplicacionCurvaEliptica(const Punto f1, const mpz_class f2, const mpz_class a, 
						const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);
						
mpz_class maxKPotenciaSuave(const mpz_class k);


mpz_class sqrt_mod(const mpz_class &a, const mpz_class &n, bool debug = false);

vector<vector<bool>> transpose(vector<vector<bool>>& matrix);

vector<vector<bool>> gaussian_elimination(const vector<vector<bool>> &matrix);

vector<vector<bool>> find_solutions(const vector<vector<bool>>& matrix);

vector<mpz_class> cocientes_fraccion_continua(const mpz_class &num, const mpz_class &den);

bool resuelveEcuacionCuadratica(mpz_class &p, mpz_class &q, mpz_class a, mpz_class b, mpz_class c);

#endif

