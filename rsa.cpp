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
bool strongPseudoprime(const mpz_t p, const mpz_t base, bool debug = false){
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
		//cout << "Estado generado: " << a << endl;
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
    //Si obtenemos un numero negativo, lo multiplicamos por -1
    if(mpz_sgn(p) <= 0){
    	mpz_mul_si(p, p, -1);
    }
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
    cout << "Buscando s,t..." << endl;
    //Generamos dos primos grandes s, t
    do{
		generate_prime(s, bits, state);
		generate_prime(t, bits, state);
		cout << "s = " << mpz_get_str (nullptr, 10, s) << endl;
		cout << "t = " << mpz_get_str (nullptr, 10, t) << endl;
		
		//Elegimos un numero aleatorio i
		mpz_urandomb(i, state, bits);
		cout << "i = " << mpz_get_str (nullptr, 10, i) << endl;
		//Calculamos r: r sea primo
		do{
			mpz_add_ui(i, i, 1);	//i++
			mpz_mul_ui(r, i, 2);	//r = 2*i*t + 1
			mpz_mul(r, r, t);
			mpz_add_ui(r, r, 1);
		}while(!millerRabin(r,10,state));
		cout << "r = " << mpz_get_str (nullptr, 10, r) << endl;
		
		//Calculamos p_0
		//p_0 = 2*s^{r-2} (mod r) *s - 1
		mpz_sub_ui(p_0, r, 2);		//p_0 = r - 2
		cout << "p_01 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_powm_sec(p_0, s, p_0, r);	//p_0 = s^(p_0) (mod r)
		cout << "p_02 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_mul_ui(p_0, p_0, 2);	//p_0 = 2*p_0
		cout << "p_03 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_fdiv_r(p_0, p_0, r);	//p_0 = p_0 (mod r)
		cout << "p_04 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_mul(p_0, p_0, s);		//p_0 = p_0*s
		cout << "p_05 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_sub_ui(p_0, p_0, 1);	//p_0 = p_0 - 1
		cout << "p_06 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		if(!mpz_odd_p(p_0))
			cout << "p_0 no es impar: repetimos el proceso" << endl;
	}while(!mpz_odd_p(p_0));	//p_0 debe ser impar, si no, repetimos el proceso
	
    //Elegimos un numero aleatorio j
    mpz_urandomb(j, state, bits);
    cout << "j = " << mpz_get_str (nullptr, 10, j) << endl;
    //aux = 2rs
    mpz_mul_ui(aux, r, 2);
    mpz_mul(aux, aux, s);
    //p = 2rsj + p_0
    mpz_mul(p, aux, j);
    cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
    cout << "p_0 = " << mpz_get_str (nullptr, 10, p_0) << endl;
    mpz_add(p, p, p_0);
    cout << "aux = " << mpz_get_str (nullptr, 10, aux) << endl;
    cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
    //Mientras p no sea primo, le sumamos aux
    while(!millerRabin(p,10,state)){
    	mpz_add(p, p, aux);
    	//cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
    }
    
    mpz_clears(p_0, r, s, t, i, j, nullptr);
}

// Función principal para generar claves RSA
void generate_rsa_key(int bits, mpz_t n, mpz_t e, mpz_t d, gmp_randstate_t state, bool strong_prime = true) {
	mpz_t p, q, phi_n, mcd;
	mpz_inits(p, q, phi_n, mcd, nullptr);
	unsigned int e_default = 65537;
	//Elegimos dos primos p y q aleatoriamente
	if(strong_prime){
		cout << "Buscando p..." << endl;
		generate_strong_prime(p, bits, state);
		cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
		cout << "Buscando q..." << endl;
		generate_strong_prime(q, bits, state);
		cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
		
	}else{
		cout << "Buscando p..." << endl;
		generate_prime(p, bits, state);
		cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
		cout << "Buscando q..." << endl;
		generate_prime(q, bits, state);
		cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
	}
	//Calculamos n y phi(n)
	cout << "Calculando n y phi(n)..." << endl;
	mpz_mul(n, p, q);		//n = p*q
	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_mul(phi_n,p, q);
	mpz_add_ui(p, p, 1);
	mpz_add_ui(q, q, 1);
	cout << "n = " << mpz_get_str (nullptr, 10, n) << endl;
	cout << "phi(n) = " << mpz_get_str (nullptr, 10, phi_n) << endl;
	if(mpz_sgn(phi_n) <= 0){
		exit(1);
	}
	cout << "Calculando e..." << endl;
	//Elegimos el exponente de cifrado e = 2^16-1 = 65.535 si < phi_n
	if(mpz_cmp_ui(phi_n, e_default) > 0){
		cout << "Asignando e_default..." << endl;
		mpz_set_ui(e, e_default);
	}
	//En otro caso, lo elegimos aleatoriamente (un bit mas pequeño que phi_n)
	else{
		cout << "Asignando e aleatorio..." << endl;
		mpz_urandomb(e, state, mpz_sizeinbase(phi_n, 10)-1);
		if(!mpz_odd_p(e)){
			mpz_sub_ui(e, e, 1);
		}
	}
	
	mpz_gcd(mcd, e, phi_n);
	cout << "mcd(e,phi_n) = " << mpz_get_str (nullptr, 10, mcd) << endl;
	//Mientras e y phi_n no sean coprimos
	while(mpz_cmp_ui(mcd, 1) > 0 || mpz_cmp_ui(e, 3) < 0 || mpz_cmp(e, phi_n) > 0 ){
		
		if(mpz_cmp_ui(mcd, 1) > 0)
			cout << "Asignando e aleatorio (eran coprimos)..." << endl;
		if(mpz_cmp_ui(e, 3) < 0)
			cout << "Asignando e aleatorio (e era menor que 3)..." << endl;
		if(mpz_cmp(e, phi_n) > 0){
			cout << "Asignando e aleatorio (e era mayor que phi(n))..." << endl;
			cout << "phi(n) = " << mpz_get_str (nullptr, 10, phi_n) << endl;
		}
		mpz_urandomb(e, state, mpz_sizeinbase(phi_n, 10)-1);
		if(!mpz_odd_p(e)){
			mpz_sub_ui(e, e, 1);
		}
		
		mpz_gcd(mcd, e, phi_n);
		cout << "e      = " << mpz_get_str (nullptr, 10, e) << endl;
		cout << "mcd(e,phi_n) = " << mpz_get_str (nullptr, 10, mcd) << endl;
	}
	
	cout << "Calculando d..." << endl;
	
	//Calculamos el inverso de e en Z_phi(n) = d
	mpz_invert(d, e, phi_n);
	
	//TODO: Comprobar ataque Wiener
	
	mpz_clears(p, q, phi_n, mcd, nullptr);
}

void cifra_RSA(mpz_t m_cifrado, mpz_t m, mpz_t e, mpz_t n){
	//Hacemos la operacion m_cifrado = m^e (mod n)
	mpz_powm(m_cifrado, m, e, n);
}

int main() {
	//Inicializacion de numeros aleatorios
	gmp_randstate_t state;
	gmp_randinit_default(state);
    gmp_randseed_ui(state, static_cast<unsigned long>(time(nullptr)));
    
    //Programa principal
	bool pseudoprimo;
	long unsigned int r;
	/*
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
	*/
	//Prueba RSA
	mpz_t n, e, d, prueba;
	mpz_inits(n, e, d, prueba, nullptr);
	
	//Generamos las claves de RSA
	generate_rsa_key(8, n, e, d, state, true);
	
	//Generamos un mensaje de prueba que enviar
	mpz_set_ui(prueba, 5);
	
	cout << "n = " << mpz_get_str (nullptr, 10, n) << endl;
	cout << "e = " << mpz_get_str (nullptr, 10, e) << endl;
	cout << "d = " << mpz_get_str (nullptr, 10, d) << endl;
	cout << "prueba = " << mpz_get_str (nullptr, 10, prueba) << endl;
	
	cifra_RSA(prueba, prueba, e, n);
	cout << "criptograma = " << mpz_get_str (nullptr, 10, prueba) << endl;
	cifra_RSA(prueba, prueba, d, n);
	
	cout << "prueba = m^(e*d) (mod n) = " << mpz_get_str (nullptr, 10, prueba) << endl;
	
	mpz_clears(n, e, d, prueba, nullptr);
	
    return 0;
}

