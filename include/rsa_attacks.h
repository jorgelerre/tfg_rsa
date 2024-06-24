#ifndef RSAATTACKS_H
#define RSAATTACKS_H

#include <gmp.h>
#include <gmpxx.h>

using namespace std;

/**
 * @file RSAAttacks.h
 * @brief Contiene funciones para realizar diversos ataques criptográficos contra RSA.
 */

/**
 * Aplica el ataque de Fermat para factorizar un número n.
 * @param n Número entero a factorizar.
 * @param state Estado del generador de números aleatorios GMP.
 * @return El factor encontrado de n.
 */
mpz_class ataqueFermat(mpz_class n, gmp_randstate_t state);

/**
 * Aplica el ataque de Kraitchik para factorizar un número n.
 * @param n Número entero a factorizar.
 * @param state Estado del generador de números aleatorios GMP.
 * @return El factor encontrado de n.
 */
mpz_class ataqueKraitchik(mpz_class n, gmp_randstate_t state);

/**
 * Aplica el algoritmo Rho de Pollard para factorizar un número n.
 * @param n Número entero a factorizar.
 * @param state Estado del generador de números aleatorios GMP.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return El factor encontrado de n.
 */
mpz_class rhoPollard(mpz_class n, gmp_randstate_t state, bool debug = false);

/**
 * Aplica el algoritmo p-1 de Pollard para factorizar un número n.
 * @param n Número entero a factorizar.
 * @param state Estado del generador de números aleatorios GMP.
 * @param k Parámetro k para el algoritmo.
 * @param att Número de intentos antes de detenerse.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return El factor encontrado de n.
 */
mpz_class p1Pollard(mpz_class n, gmp_randstate_t state, const mpz_class k, 
                    const unsigned int att, bool debug = true);

/**
 * Realiza la factorización de n utilizando curvas elípticas.
 * @param n Número entero a factorizar.
 * @param state Estado del generador de números aleatorios GMP.
 * @param k Parámetro k para el algoritmo.
 * @param att Número máximo de iteraciones para el algoritmo.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return El factor encontrado de n.
 */
//mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k, 
//                                       const unsigned int att = 1000, bool debug = false);

/**
 * Realiza la factorización de n utilizando la criba cuadrática.
 * @param n Número entero a factorizar.
 * @param k Parámetro k para el algoritmo.
 * @param tam_tabla Tamaño de la tabla para la criba.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return El factor encontrado de n.
 */
mpz_class factorizacionCribaCuadratica(const mpz_class &n, const mpz_class &k, 
                                       const unsigned int tam_tabla, bool debug = false);

/**
 * Realiza el ataque de Wiener para encontrar la clave privada d dado e y n.
 * @param d Referencia a la clave privada a encontrar.
 * @param p Referencia a un factor primo de n.
 * @param q Referencia al otro factor primo de n.
 * @param e Clave pública.
 * @param n Módulo de la clave RSA.
 * @param debug Habilita la salida de depuración si es verdadero.
 * @return Verdadero si el ataque tiene éxito, falso si falla.
 */
bool ataqueWiener(mpz_class &d, mpz_class &p, mpz_class &q, 
                  const mpz_class e, const mpz_class n, bool debug = false);

#endif // RSAATTACKS_H

