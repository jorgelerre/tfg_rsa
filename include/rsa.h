#ifndef RSA_H
#define RSA_H

#include <gmp.h>
#include <gmpxx.h>

/**
 * @file RSA.h
 * @brief Contiene funciones para la generación de claves RSA y operaciones de cifrado y descifrado RSA.
 */

/**
 * Genera una clave RSA pública y privada.
 * @param n Referencia al módulo RSA generado.
 * @param e Referencia al exponente público generado.
 * @param d Referencia al exponente privado generado.
 * @param bits Número de bits deseado para el módulo n.
 * @param state Estado del generador de números aleatorios GMP.
 * @param strong_prime Indica si se debe generar un primo robusto para p y q.
 * @param low_d Indica si se desea un valor bajo para el exponente d (opcional).
 * @param debug Habilita la salida de depuración si es verdadero (opcional).
 */
void generaClaveRSA(mpz_class &n, mpz_class &e, mpz_class &d, unsigned int bits, 
                      gmp_randstate_t state, bool strong_prime, bool low_d = false, bool debug = false);

/**
 * Realiza el cifrado RSA de un mensaje m.
 * @param m Mensaje a cifrar.
 * @param e Exponente de cifrado.
 * @param n Módulo RSA.
 * @return El mensaje cifrado.
 */
mpz_class cifraRSA(mpz_class m, mpz_class e, mpz_class n);

/**
 * Realiza el descifrado RSA de un mensaje cifrado m.
 * @param m Mensaje cifrado a descifrar.
 * @param d Exponente de descifrado.
 * @param n Módulo RSA.
 * @return El mensaje descifrado.
 */
mpz_class descifraRSA(mpz_class m, mpz_class d, mpz_class n);

#endif // RSA_H

