#include "utils.h"
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;


/*
 * Genera un numero aleatorio corto
 * @it Numero de iteraciones que aplicar la prueba
 * @return Un numero aleatorio de 64 bits
*/
unsigned long int generate_random_limb(gmp_randstate_t state) {
    mpz_t random_num;
    unsigned long int r;
    mpz_init(random_num);

    // Genera un número aleatorio de tamaño máximo de un unsigned long int
    mpz_urandomb(random_num, state, GMP_LIMB_BITS);

    // Convierte el número generado a unsigned long int
    r = mpz_getlimbn(random_num, 0);

    // Limpieza
    mpz_clear(random_num);
    
    return r;
}


/*
 * Comprueba si un numero es pseudoprimo fuerte en cierta base
 * @p Numero entero a analizar
 * @base Base donde estudiar la pseudoprimalidad de p
 * @return True si p es pseudoprimo fuerte en base, False en otro caso.
 */
bool strongPseudoprime(const mpz_class p, const mpz_class base, bool debug){
	bool is_pseudoprime = false;
	mpz_class mcd;
	mpz_class q, r;
	mpz_class t, s;
	mpz_class b;
	
	//Comprobamos que sean coprimos y que el primo que comprobamos sea mayor que 1
	mpz_gcd(mcd.get_mpz_t(), p.get_mpz_t(), base.get_mpz_t());
	
	if(debug)
		cout << "mcd de a y la base = " << mcd << endl;
		
	if(mcd == 1 && p > 1){
		//Calculamos s,t tales que n-1 = 2^s*t con t impar
		t = p - 1;			// t = p - 1
		s = 0;				// s = 0
		mpz_fdiv_qr_ui(q.get_mpz_t(), r.get_mpz_t(), t.get_mpz_t(), 2);		// Dividimos n-1 entre 2, obteniendo cociente q y resto r
		while(r == 0){	// Si el resto es 0, t es impar aun
			t = q;				// t = q
			s = s + 1;		// s = s + 1
			mpz_fdiv_qr_ui(q.get_mpz_t(), r.get_mpz_t(), t.get_mpz_t(), 2);	// (q,r) = t / 2, t % 2
		}
		if(debug){
			cout << "s = " << s << endl;
			cout << "t = " << t << endl;
		}
		
		//Una vez obtenidos s,t, calculamos la primera condicion: (base^t (mod p) == 1)
		mpz_powm(b.get_mpz_t(), base.get_mpz_t(), t.get_mpz_t(), p.get_mpz_t());	// b = base^t mod p
		q = p - 1;		// q = p - 1
		if(b == 1 || b == q){	// Si b = 1 o b = p -1 --> PSEUDOPRIMO
			is_pseudoprime = true;
		}
		
		//Calculamos la segunda condicion: base^{2^e*r} % p = p - 1
		r = 1;	// r = 1
		//Mientras s > r, b 
		while(s > r && !is_pseudoprime){
			b = (b*b) % p;
			if(b == p - 1){	
				is_pseudoprime = true;
			}
			r++;
		}
	}
	
	return is_pseudoprime;
}

/*
 * Aplica el test de Miller Rabin sobre un primo cierto numero de veces
 * @p Numero entero a analizar
 * @it Numero de iteraciones que aplicar la prueba
 * @state El estado del generador de numeros aleatorios
 * @return True si p es un primo con confianza 1-4^{-it}, false en otro caso.
*/
bool millerRabin(const mpz_class p, int it, gmp_randstate_t state){
	bool no_fail = true;
	mpz_class b;
	
	//Realizamos el test de las divisiones sucesivas con los primeros primos (<8 bits)
	const int NUM_PRIMES = 54;
	const unsigned int primes[NUM_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
			   59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
	 		   157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};
	
	for(int i = 0; i < NUM_PRIMES && no_fail; i++){
		b = p % primes[i];
		if(b == 0)
			no_fail = false;
	}
	
	//Realizamos it veces la comprobacion de Miller-Rabin
	for(int i = 0; i < it && no_fail; i++){
		//Generamos un numero aleatorio (pequeño (1 miembro) para aligerar los calculos)
		b = generate_random_limb(state);
		//Comprobamos la pseudoprimalidad en dicha base
		no_fail = strongPseudoprime(p, b);
	}
	
	return no_fail;
}

/**
 * Función para generar un número primo de 'bits' bits
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 * @return Un primo del tamanio especificado
 */
mpz_class generate_prime(unsigned int bits, gmp_randstate_t state) {
    mpz_class p;
    
    do{
		mpz_rrandomb(p.get_mpz_t(), state, bits); // Genera un número aleatorio de 'bits' bits
		//Si obtenemos un numero negativo, lo multiplicamos por -1
		if(p < 0){
			p *= -1;
		}
		//Si el numero es par, le sumamos 1
		if(p % 2 == 0){
			p++;
		}
		//Mientras p no pase el test de Miller-Rabin, avanzamos al siguiente numero impar
		while(!millerRabin(p, 10, state)){
			p += 2;
		}
    }while(mpz_sizeinbase(p.get_mpz_t(), 2) != bits);	//Comprobamos que p no se pase del tamanio especificado
    
    return p;
}

/**
 * Función para generar un número primo robusto de 'bits' bits
 * 
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 * @return Dato de tipo mpz_class con el primo generado
 */
mpz_class generate_strong_prime(unsigned int bits, gmp_randstate_t state, bool debug) {
    mpz_class p,p_0,r,s,t,i,j,aux;
    unsigned int bits_j, bits_1, bits_2;
    
    if(debug) cout << "Buscando s,t..." << endl;
    //Generamos dos primos grandes s, t
    bits_1 = (bits-log2(bits))/2 - 4;
    bits_2 = bits_1 - log2(bits_1) - 7;
    if(bits_1 < 4 || bits_2 < 4){
    	bits_1 = bits / 2;
    	bits_2 = bits / 2 - 1;
    }
    do{
		s = generate_prime(bits_1, state);
		t = generate_prime(bits_2, state);
		if(debug){
			cout << "s = " << s << endl;
			cout << "t = " << t << endl;
		}
		//Elegimos un numero aleatorio i
		mpz_urandomb(i.get_mpz_t(), state, log2(bits));
		if(debug) cout << "i = " << i << endl;
		//Calculamos r: r sea primo
		do{
			i++;
			r = 2*i*t + 1;
		}while(!millerRabin(r,10,state));
		if(debug){	
			cout << "r = " << r << endl;
			cout << "tamanio r = " << mpz_sizeinbase(r.get_mpz_t(), 2) << endl;
		}
		//Calculamos p_0
		//p_0 = ((2*s^{r-2}) % r) *s - 1
		p_0 = r - 2;
		mpz_powm_sec(p_0.get_mpz_t(), s.get_mpz_t(), 
					 p_0.get_mpz_t(), r.get_mpz_t());	//p_0 = s^(p_0) (mod r)
		//cout << "p_02 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		p_0 = ((2*p_0) % r)*s - 1;

		//if(!mpz_odd_p(p_0))
			//cout << "p_0 no es impar: repetimos el proceso" << endl;
	}while(p_0 % 2 == 0);	//p_0 debe ser impar, si no, repetimos el proceso
	if(debug){
		cout << "tamanio s = " << mpz_sizeinbase(s.get_mpz_t(), 2) << endl;
		cout << "tamanio t = " << mpz_sizeinbase(t.get_mpz_t(), 2) << endl;
		cout << "tamanio p_0 = " << mpz_sizeinbase(p_0.get_mpz_t(), 2) << endl;
    }
    //aux = 2rs
    aux = 2*r*s;
    
    //Elegimos un numero aleatorio j tal que tenga los bits suficientes para que
    //p tenga 'bits' bits.
    if(debug){
		cout << "Bits aux: " << mpz_sizeinbase(aux.get_mpz_t(), 2) << endl;
		cout << "Bits: " << bits << endl;
	}
    bits_j = bits - mpz_sizeinbase(aux.get_mpz_t(), 2);
    if(bits_j < 1){
    	bits_j = 1;
    }
    mpz_urandomb(j.get_mpz_t(), state, bits_j);
    if(debug) cout << "j = " << j << endl;
    
    //p = 2rsj + p_0
    p = aux*j + p_0; 
    if(debug) {
		cout << "aux = " << aux << endl;
		cout << "j = " << j << endl;
		cout << "p = " << p << endl;
    }
    //Mientras p no sea primo, le sumamos aux
    while(!millerRabin(p,10,state)){
    	p += aux;
    }
    if(debug) cout << "p final = " << p << endl;
    return p;
}

