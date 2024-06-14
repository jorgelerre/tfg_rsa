#include "rsa.h"
#include "rsa_attacks.h"
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

int main() {
	//Inicializacion de numeros aleatorios
	gmp_randstate_t state;
	gmp_randinit_default(state);
    gmp_randseed_ui(state, static_cast<unsigned long>(time(nullptr)));
    
    //Programa principal - Menu
	
	
	/*
	bool pseudoprimo;
	long unsigned int r;
	mpz_class p, base;
	mpz_inits(p, base, nullptr);
	mpz_set_ui(p, 503);
	mpz_set_ui(base, 2);
	
	pseudoprimo = millerRabin(p, 5, state);
	cout << "Pseudoprimo: " << pseudoprimo << endl;
	
	generate_prime(p, 32, state);
	r = mpz_getlimbn(p, 0);
	cout << "Primo encontrado: " << r << endl;
	
	mpz_clears(p, base, state, nullptr);
	
	//Prueba RSA
	mpz_t n, e, d, prueba;
	mpz_inits(n, e, d, prueba, nullptr);
	
	//Generamos las claves de RSA
	generate_rsa_key(30, n, e, d, state, true);
	
	//Generamos un mensaje de prueba que enviar
	mpz_set_ui(prueba, 5);
	
	cout << "n = " << mpz_get_str (nullptr, 10, n) << endl;
	cout << "bits n = " << mpz_sizeinbase(n, 2) << endl;
	cout << "e = " << mpz_get_str (nullptr, 10, e) << endl;
	cout << "d = " << mpz_get_str (nullptr, 10, d) << endl;
	cout << "prueba = " << mpz_get_str (nullptr, 10, prueba) << endl;
	
	cifra_RSA(prueba, prueba, e, n);
	cout << "criptograma = " << mpz_get_str (nullptr, 10, prueba) << endl;
	cifra_RSA(prueba, prueba, d, n);
	
	cout << "prueba = m^(e*d) (mod n) = " << mpz_get_str (nullptr, 10, prueba) << endl;
	
	mpz_clears(n, e, d, prueba, nullptr);
	*/
	// Prueba suma en curvas elipticas
	/*
	mpz_t a, b, p, inv;
	Punto suma, s1, s2;
	int res = 0;
	mpz_inits(suma.x, suma.y, s1.x, s1.y, s2.x, s2.y, a, b, p, inv, nullptr);
	
	mpz_set_ui(s1.x, 9696);
	mpz_set_ui(s1.y, 506);
	
	mpz_set_ui(s2.x, 7878);
	mpz_set_ui(s2.y, 10200);
	
	mpz_set_ui(a, 1);
	mpz_set_ui(b, 1);
	mpz_set_ui(p, 10403);
	
	
	res = sumaCurvaEliptica(suma, s1, s2, a, b, p, inv);
	cout << "res = " << res << endl;
	cout << "suma = [" << mpz_get_str (nullptr, 10, suma.x) << ", " 
		 		       << mpz_get_str (nullptr, 10, suma.y) << "]" << endl;
	cout << "s1 = [" << mpz_get_str (nullptr, 10, s1.x) << ", " 
		 		       << mpz_get_str (nullptr, 10, s1.y) << "]" << endl;
	cout << "s2 = [" << mpz_get_str (nullptr, 10, s2.x) << ", " 
		 		       << mpz_get_str (nullptr, 10, s2.y) << "]" << endl;
	cout << "a = " << mpz_get_str (nullptr, 10, a) << endl;
	cout << "b = " << mpz_get_str (nullptr, 10, b) << endl;
	cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
	cout << "inv = " << mpz_get_str (nullptr, 10, inv) << endl;
	
	mpz_clears(suma.x, suma.y, s1.x, s1.y, s2.x, s2.y, a, b, p, inv, nullptr);
	*/
	// Prueba multiplicacion en curvas elipticas
	/*
	mpz_t a, b, p, inv, f2, k;
	Punto mul, f1;
	mpz_inits(mul.x, mul.y, f1.x, f1.y, f2, a, b, p, inv, k, nullptr);
	mpz_set_ui(f1.x, 0);
	mpz_set_ui(f1.y, 1);
	mpz_set_ui(k,20);
	maxKPotenciaSuave(f2, k);
	//mpz_set_ui(f2, 1);
	mpz_set_ui(a, 1);
	mpz_set_ui(b, 1);
	mpz_set_ui(p, 10403);
	
	multiplicacionCurvaEliptica(mul, f1, f2, a, b, p, inv);
	cout << "mul = [" << mpz_get_str (nullptr, 10, mul.x) << ", " 
					  << mpz_get_str (nullptr, 10, mul.y) << "]" << endl;
	cout << "f1 = [" << mpz_get_str (nullptr, 10, f1.x) << ", " 
				     << mpz_get_str (nullptr, 10, f1.y) << "]" << endl;
	cout << "f2 = " << mpz_get_str (nullptr, 10, f2) << endl;
	cout << "a = " << mpz_get_str (nullptr, 10, a) << endl;
	cout << "b = " << mpz_get_str (nullptr, 10, b) << endl;
	cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
	cout << "inv = " << mpz_get_str (nullptr, 10, inv) << endl;
	*/
	// Prueba Curvas Elipticas
	unsigned int att = 1000;
	mpz_class n, p, q, k;
	
	//mpz_set_ui(n, 265905204186593);
	n = 265905204186593;
	k = 30;
	p = factorizacionCurvasElipticas(n, state, k, att);
	if(p != n){
		cout << "n = " << n << endl;
		cout << "p = " << p << endl;
		cout << "q = " << n / p << endl;
	}
	else{
		cout << "No se ha conseguido factorizar el numero :(" << endl;
	}
	
	mpz_class sqrt_a;
	
	mpz_class n1, p1, q1, k1;
	n1 = 265905204186593;
	//n1 = 15167;
	k1 = 1000;
	
	//sqrt_mod(sqrt_a, k, n);
	//cout << sqrt_a << endl;
	p = factorizacionCribaCuadratica(n1, k1, 1000000);
	
    return 0;
}

