#ifndef RSA_H
#define RSA_H

#include <gmp.h>
#include <gmpxx.h>

void generate_rsa_key(mpz_class &n, mpz_class &e, mpz_class &d, unsigned int bits, 
				      gmp_randstate_t state, bool strong_prime, bool low_d = false, bool debug = false);
mpz_class cifra_RSA(mpz_class m, mpz_class e, mpz_class n);
mpz_class descifra_RSA(mpz_class m, mpz_class d, mpz_class n);
#endif

