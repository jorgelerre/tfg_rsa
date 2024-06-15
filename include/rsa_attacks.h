#ifndef RSAATTACKS_H
#define RSAATTACKS_H

#include <gmp.h>
#include <gmpxx.h>

using namespace std;

mpz_class ataqueFermat(mpz_class n, gmp_randstate_t state);

mpz_class ataqueKraitchik(mpz_class n, gmp_randstate_t state);

mpz_class rhoPollard(mpz_class n, gmp_randstate_t state, bool debug = false);

mpz_class p1Pollard(mpz_class n, gmp_randstate_t state, const mpz_class k, 
					const unsigned int att, bool debug = true);

mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k, 
									   const unsigned int att = 1000, bool debug = false);

mpz_class factorizacionCribaCuadratica(const mpz_class &n, const mpz_class &k, 
									   const unsigned int tam_tabla, bool debug = false);

bool ataqueWiener(mpz_class &d, mpz_class &p, mpz_class &q, 
				  const mpz_class e, const mpz_class n, bool debug = false);

#endif

