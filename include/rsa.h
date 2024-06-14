#ifndef RSA_H
#define RSA_H

#include <gmp.h>
#include <gmpxx.h>

void generate_rsa_key(unsigned bits, mpz_class &n, mpz_class &e, mpz_class &d, gmp_randstate_t state, 
				      bool strong_prime = true, bool debug = false);
mpz_class cifra_RSA(mpz_class m, mpz_class e, mpz_class n);
mpz_class descifra_RSA(mpz_class m, mpz_class d, mpz_class n);
#endif

