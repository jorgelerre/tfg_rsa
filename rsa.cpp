#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>

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
bool strongPseudoprime(const mpz_t p, const mpz_t base, bool debug = true){
	cout << "HOLA" << endl << flush;
	bool is_pseudoprime = false;
	mpz_t mcd;
	mpz_t q, r;
	mpz_t t, s;
	mpz_t b;
	mpz_inits(mcd, q, r, t, s, b, nullptr);
	
	//Comprobamos que sean coprimos y que el primo que comprobamos sea mayor que 1
	mpz_gcd(mcd, p, base);
	if(debug)
		cout << "mcd de a y la base = " << mpz_get_str (nullptr, 10, mcd) << endl;
	if(mpz_cmp_ui(mcd,1) == 0 && mpz_cmp_ui(p,1) > 0){
		//Calculamos s,t tales que n-1 = 2^s*t con t impar
		mpz_sub_ui(t, p, 1);			// t = p - 1
		mpz_set_ui(s, 0);				// s = 0
		mpz_fdiv_qr_ui(q, r, t, 2);		// Dividimos n-1 entre 2, obteniendo cociente q y resto r
		while(mpz_cmp_ui(r,0) == 0){	// Si el resto es 0, t es impar aun
			mpz_set(t, q);				// t = q
			mpz_add_ui(s, s, 1);		// s = s + 1
			mpz_fdiv_qr_ui(q, r, t, 2);	// (q,r) = t / 2, t % 2
		}
		if(debug){
			cout << "s = " << s << endl;
			cout << "t = " << t << endl;
		}
		//Una vez obtenidos s,t, calculamos la primera condicion: (base^t (mod p) == 1)
		mpz_powm(b, base, t, p);	// b = base^t mod p
		mpz_sub_ui(q, p, 1);		// q = p - 1
		if(mpz_cmp_ui(b,1) == 0 || mpz_cmp(b,q) == 0){	// Si b = 1 o b = p -1 --> PSEUDOPRIMO
			is_pseudoprime = true;
		}
		
		//Calculamos la segunda condicion: base^{2^e*r} % p = p - 1
		mpz_set_ui(r, 1);	// r = 1
		//Mientras s > r, b 
		while(mpz_cmp(s,r) > 0 && !is_pseudoprime){
			// b = (b*b) % p;
			mpz_mul(b,b,b);
			mpz_fdiv_r(b,b,p);
			
			if(mpz_cmp(b,q) == 0){	// b == p - 1
				is_pseudoprime = true;
			}
			mpz_add_ui(r, r, 1);	//r++
		}
	}
	
	mpz_clears(mcd, q, r, t, s, b, nullptr);
	return is_pseudoprime;
}

/*
 * Aplica el test de Miller Rabin sobre un primo cierto numero de veces
 * @p Numero entero a analizar
 * @it Numero de iteraciones que aplicar la prueba
 * @state El estado del generador de numeros aleatorios
 * @return True si p es un primo con confianza 1-4^{-it}, false en otro caso.
*/
bool millerRabin(mpz_t p, int it, gmp_randstate_t state){
	bool no_fail = true;
	unsigned long int a;
	mpz_t b;
	
	mpz_init(b);
	
	//Realizamos it veces la comprobacion de Miller-Rabin
	for(int i = 0; i < it && no_fail; i++){
		//Generamos un numero aleatorio (pequeño (1 miembro) para aligerar los calculos)
		a = generate_random_limb(state);
		cout << "Estado generado: " << a << endl;
		mpz_set_ui(b,a);
		//Comprobamos la pseudoprimalidad en dicha base
		no_fail = strongPseudoprime(p, b);
	}
	
	mpz_clear(b);
	
	return no_fail;
}


// Función para generar un número primo de 'bits' bits
void generate_prime(mpz_t p, int bits, gmp_randstate_t state) {
    mpz_rrandomb(p, state, bits); // Genera un número aleatorio de 'bits' bits
    //Si el numero es par, le sumamos 1
    if(!mpz_odd_p(p)){
    	mpz_add_ui(p, p, 1);
    }
    //Mientras p no pase el test de Miller-Rabin, avanzamos al siguiente numero impar
    while(!millerRabin(p, 10, state)){
    	mpz_add_ui(p, p, 2);
    }
}

// Función para generar un número primo robusto de 'bits' bits
void generate_strong_prime(mpz_t p, int bits, gmp_randstate_t state) {
    mpz_t p_0,r,s,t,i,j,aux;
    
    mpz_inits(p_0, r, s, t, i,j,aux, nullptr);
    
    //Generamos dos primos grandes s, t
    generate_prime(s, bits, state);
    generate_prime(t, bits, state);
    
    //Elegimos un numero aleatorio i
    mpz_urandomb(i, state, bits);
    
    //Calculamos r: r sea primo
    do{
    	mpz_add_ui(i, i, 1);	//i++
    	mpz_mul_ui(r, i, 2);	//r = 2*i*t + 1
    	mpz_mul(r, r, t);
    	mpz_add_ui(r, r, 1);
    }while(!millerRabin(r,10,state));
    
    //Calculamos p_0
    //p_0 = 2*s^{r-2} (mod r) *s - 1
    mpz_sub_ui(p_0, r, 2);
    mpz_powm(p_0, s, p_0, r);
    mpz_mul_ui(p_0, p_0, 2);
    mpz_fdiv_r(p_0, p_0, r);
    mpz_mul(p_0, p_0, s);
    mpz_sub_ui(p_0, p_0, 1);
    
    //Elegimos un numero aleatorio j
    mpz_urandomb(j, state, bits);
    //aux = 2rs
    mpz_mul_ui(aux, r, 2);
    mpz_mul(aux, aux, s);
    //p = 2rsj + p_0
    mpz_mul(p, aux, j);
    mpz_add(p, p, p_0);
    
    //Mientras p no sea primo, le sumamos aux
    while(!millerRabin(p,10,state)){
    	mpz_add(p, p, aux);
    }
    
    mpz_clears(p_0, r, s, t, i, j, nullptr);
}

// Función principal para generar claves RSA
// TODO: Implementar parametro para elegir generacion de primos (robustos o no)
void generate_rsa_key(int bits, mpz_t n, mpz_t e, mpz_t d, gmp_randstate_t state) {
	mpz_t p, q, phi_n, mcd;
	mpz_inits(p, q, phi_n, mcd, nullptr);
	unsigned int e_default = 65535;
	//Elegimos dos primos p y q aleatoriamente
	generate_prime(p, bits, state);
	generate_prime(q, bits, state);
	
	//Calculamos n y phi(n)
	mpz_mul(n,p,q);		//n = p*q
	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_mul(phi_n,p, q);
	mpz_add_ui(p, p, 1);
	mpz_add_ui(q, q, 1);
	
	//Elegimos el exponente de cifrado e = 2^16-1 = 65.535 si < phi_n
	if(mpz_cmp_ui(phi_n, e_default) > 0){
		mpz_set_ui(e, e_default);
	}
	//En otro caso, lo elegimos aleatoriamente (un bit mas pequeño que phi_n)
	else{
		mpz_urandomb(e, state, mpz_sizeinbase(phi_n, 2));
		if(!mpz_odd_p(e)){
			mpz_sub_ui(e, e, 1);
		}
	}
	
	mpz_gcd(mcd, e, phi_n);
	mpz_sub_ui(mcd, mcd, 1);
	
	//Mientras e y phi_n no sean coprimos
	while(mpz_cmp_ui(mcd, 0) > 0 && mpz_cmp_ui(e, 2) > 0 && mpz_cmp_ui(e, 2) > 0 ){
		mpz_urandomb(e, state, mpz_sizeinbase(phi_n, 2));
		if(!mpz_odd_p(e)){
			mpz_sub_ui(e, e, 1);
		}
		
		mpz_gcd(mcd, e, phi_n);
		mpz_sub_ui(mcd, mcd, 1);
	}
	
	//Calculamos el inverso de e en Z_phi(n) = d
	mpz_invert(d, e, phi_n);
	
	//TODO: Comprobar ataque Wiener
	
	mpz_clears(p, q, phi_n, mcd, nullptr);
}

int main() {
	//Inicializacion de numeros aleatorios
	gmp_randstate_t state;
	gmp_randinit_default(state);
    gmp_randseed_ui(state, static_cast<unsigned long>(time(nullptr)));
    
    //Programa principal
	bool pseudoprimo;
	long unsigned int r;
	mpz_t p, base;
	mpz_inits(p, base, nullptr);
	mpz_set_ui(p, 503);
	mpz_set_ui(base, 2);
	pseudoprimo = millerRabin(p, 5, state);
	cout << "Pseudoprimo: " << pseudoprimo << endl;
	
	
	generate_prime(p, 32, state);
	r = mpz_getlimbn(p, 0);
	cout << "Primo encontrado: " << r << endl;
	
	mpz_clears(p, base, state, nullptr);
    return 0;
}

