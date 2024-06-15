#include "rsa.h"
#include "utils.h"
#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;


// Función principal para generar claves RSA
void generate_rsa_key(mpz_class &n, mpz_class &e, mpz_class &d, unsigned int bits, gmp_randstate_t state, bool strong_prime, bool low_d, bool debug) {
    mpz_class p, q, phi_n, mcd;
    
    unsigned int e_default = 65537;

    // Elegimos dos primos p y q aleatoriamente
    if (strong_prime) {
        if (debug) cout << "Buscando p..." << endl;
        p = generate_strong_prime(bits / 2, state);
        if (debug) {
            cout << "p = " << p << endl;
            cout << "tamaño p = " << mpz_sizeinbase(p.get_mpz_t(), 2) << endl;
            cout << "Buscando q..." << endl;
        }
        q = generate_strong_prime(bits / 2, state);
        if (debug) {
            cout << "q = " << q << endl;
            cout << "Tamaño q = " << mpz_sizeinbase(q.get_mpz_t(), 2) << endl;
        }
    } else {
        if (debug) cout << "Buscando p..." << endl;
        p = generate_prime(bits / 2, state);
        if (debug) {
            cout << "p = " << p << endl;
            cout << "tamaño p = " << mpz_sizeinbase(p.get_mpz_t(), 2) << endl;
            cout << "Buscando q..." << endl;
        }
        q = generate_prime(bits / 2, state);
        if (debug) {
            cout << "q = " << q << endl;
            cout << "Tamaño q = " << mpz_sizeinbase(q.get_mpz_t(), 2) << endl;
        }
    }

    // Calculamos n y phi(n)
    if (debug) cout << "Calculando n y phi(n)..." << endl;
    n = p * q;
    phi_n = (p - 1) * (q - 1);

    if (debug) {
        cout << "n = " << n << endl;
        cout << "phi(n) = " << phi_n << endl;
        cout << "Calculando e..." << endl;
    }

	if(!low_d){
		// Elegimos el exponente de cifrado e = 2^16-1 = 65.537 si es < phi_n
		if (phi_n > e_default) {
		    if (debug) cout << "Asignando e_default..." << endl;
		    e = e_default;
		} else {
		    if (debug) cout << "Asignando e aleatorio..." << endl;
		    mpz_urandomb(e.get_mpz_t(), state, mpz_sizeinbase(phi_n.get_mpz_t(), 10) - 1);
		    if (e % 2 == 0) {
		        e--;
		    }
		}

		mpz_gcd(mcd.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
		if (debug) cout << "mcd(e,phi_n) = " << mcd << endl;

		// Mientras e y phi_n no sean coprimos
		while (mcd > 1 || e < 3 || e > phi_n) {
		    if (mcd > 1 && debug) cout << "Asignando e aleatorio (no eran coprimos)..." << endl;
		    if (e < 3 && debug) cout << "Asignando e aleatorio (e era menor que 3)..." << endl;
		    if (e > phi_n && debug) {
	            cout << "Asignando e aleatorio (e era mayor que phi(n))..." << endl;
	            cout << "phi(n) = " << phi_n << endl;
		    }
		    mpz_urandomb(e.get_mpz_t(), state, mpz_sizeinbase(phi_n.get_mpz_t(), 2) - 1);
		    if (e % 2 == 0) {
		        e--;
		    }

		    mpz_gcd(mcd.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
		    if (debug) {
		        cout << "e = " << e << endl;
		        cout << "mcd(e,phi_n) = " << mcd << endl;
		    }
		}

		if (debug) cout << "Calculando d..." << endl;

		// Calculamos el inverso de e en Z_phi(n) = d
		mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
	}
	//Si la opcion low_d esta activada, usaremos un d vulnerable al ataque de Wiener
	else{
		//Generamos un exponente de descifrado d pequeño
		cout << "n = " << n << endl;
		mpz_class lim = sqrt(sqrt(n))/3;
		cout << "Cota d: " << lim << endl;
		do{
			mpz_urandomm(d.get_mpz_t(), state, lim.get_mpz_t());
			cout << "d generado = " << d << endl;
			mpz_gcd(mcd.get_mpz_t(), d.get_mpz_t(), phi_n.get_mpz_t());
			cout << "mcd = " << mcd << endl;
		}while(d < 2 || mcd > 1);
		// Calculamos el inverso de d en Z_phi(n) = e
		mpz_invert(e.get_mpz_t(), d.get_mpz_t(), phi_n.get_mpz_t());
	}
}


mpz_class cifra_RSA(mpz_class m, mpz_class e, mpz_class n){
	mpz_class m_cifrado;
	//Hacemos la operacion m_cifrado = m^e (mod n)
	mpz_powm(m_cifrado.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
	return m_cifrado;
}

mpz_class descifra_RSA(mpz_class m, mpz_class d, mpz_class n){
	return cifra_RSA(m, d, n);	//Se cifra de la misma forma que se descifra
}

