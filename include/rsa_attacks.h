#ifndef RSAATTACKS_H
#define RSAATTACKS_H

#include <gmp.h>
#include <gmpxx.h>
#include <vector>

using namespace std;

struct Punto{
	mpz_class x;
	mpz_class y;
};

mpz_class ataqueFermat(mpz_class n, gmp_randstate_t state);

mpz_class ataqueKraitchik(mpz_class n, gmp_randstate_t state);

mpz_class rhoPollard(mpz_class n, gmp_randstate_t state, bool debug = false);

//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
Punto sumaCurvaEliptica(const Punto s1, const Punto s2, const mpz_class a, 
					    const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);
					    
Punto multiplicacionCurvaEliptica(const Punto f1, const mpz_class f2, const mpz_class a, 
						const mpz_class b, const mpz_class p, mpz_class &inv, bool debug = false);
						
mpz_class maxKPotenciaSuave(const mpz_class k);

mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k, const unsigned int att = 1000, bool debug = false);

mpz_class sqrt_mod(const mpz_class &a, const mpz_class &n, bool debug = false);

vector<vector<bool>> transpose(vector<vector<bool>>& matrix);

vector<vector<bool>> gaussian_elimination(const vector<vector<bool>> &matrix);

vector<vector<bool>> find_solutions(const vector<vector<bool>>& matrix);

mpz_class factorizacionCribaCuadratica(const mpz_class &n, const mpz_class &k, 
									   const unsigned int tam_tabla, bool debug = false);

#endif

