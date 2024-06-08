#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
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
bool millerRabin(const mpz_t p, int it, gmp_randstate_t state){
	bool no_fail = true;
	unsigned long int a;
	mpz_t b;
	
	mpz_init(b);
	
	//Realizamos el test de las divisiones sucesivas con los primeros primos (<8 bits)
	const int NUM_PRIMES = 54;
	const unsigned int primes[NUM_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
			   59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
	 		   157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};
	
	if(mpz_cmp_ui(p, 252) > 0){
		for(int i = 0; i < NUM_PRIMES && no_fail; i++){
			mpz_fdiv_r_ui(b,p,primes[i]);
			if(mpz_cmp_ui(b,0) == 0)
				no_fail = false;
		}
	}
	
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

/**
 * Función para generar un número primo de 'bits' bits
 * @p Dato de tipo mpz_t donde se devuelve el primo generado
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 */
void generate_prime(mpz_t p, int bits, gmp_randstate_t state) {
    do{
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
    }while(mpz_sizeinbase(p, 2) != bits);
}

/**
 * Función para generar un número primo robusto de 'bits' bits
 * @p Dato de tipo mpz_t donde se devuelve el primo generado
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 */
void generate_strong_prime(mpz_t p, int bits, gmp_randstate_t state) {
    mpz_t p_0,r,s,t,i,j,aux;
    int bits_j, bits_1, bits_2;
    mpz_inits(p_0, r, s, t, i, j, aux, nullptr);
    //cout << "Buscando s,t..." << endl;
    //Generamos dos primos grandes s, t
    bits_1 = (bits-log2(bits))/2 - 4;
    bits_2 = bits_1 - log2(bits_1) - 7;
    if(bits_1 < 4 || bits_2 < 4){
    	bits_1 = bits / 2;
    	bits_2 = bits / 2 - 1;
    }
    do{
		generate_prime(s, bits_1, state);
		generate_prime(t, bits_2, state);
		//cout << "s = " << mpz_get_str (nullptr, 10, s) << endl;
		//cout << "t = " << mpz_get_str (nullptr, 10, t) << endl;
		
		//Elegimos un numero aleatorio i
		mpz_urandomb(i, state, log2(bits));
		//cout << "i = " << mpz_get_str (nullptr, 10, i) << endl;
		//Calculamos r: r sea primo
		do{
			mpz_add_ui(i, i, 1);	//i++
			mpz_mul_ui(r, i, 2);	//r = 2*i*t + 1
			mpz_mul(r, r, t);
			mpz_add_ui(r, r, 1);
		}while(!millerRabin(r,10,state));
		//cout << "r = " << mpz_get_str (nullptr, 10, r) << endl;
		//cout << "tamanio r = " << mpz_sizeinbase(r, 2) << endl;
		//Calculamos p_0
		//p_0 = 2*s^{r-2} (mod r) *s - 1
		mpz_sub_ui(p_0, r, 2);		//p_0 = r - 2
		//cout << "p_01 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_powm_sec(p_0, s, p_0, r);	//p_0 = s^(p_0) (mod r)
		//cout << "p_02 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_mul_ui(p_0, p_0, 2);	//p_0 = 2*p_0
		//cout << "p_03 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_fdiv_r(p_0, p_0, r);	//p_0 = p_0 (mod r)
		//cout << "p_04 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_mul(p_0, p_0, s);		//p_0 = p_0*s
		//cout << "p_05 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		mpz_sub_ui(p_0, p_0, 1);	//p_0 = p_0 - 1
		//cout << "p_06 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		//if(!mpz_odd_p(p_0))
			//cout << "p_0 no es impar: repetimos el proceso" << endl;
	}while(!mpz_odd_p(p_0));	//p_0 debe ser impar, si no, repetimos el proceso
	cout << "tamanio s = " << mpz_sizeinbase(s, 2) << endl;
	cout << "tamanio t = " << mpz_sizeinbase(t, 2) << endl;
	cout << "tamanio p_0 = " << mpz_sizeinbase(p_0, 2) << endl;
    
    //aux = 2rs
    mpz_mul_ui(aux, r, 2);
    mpz_mul(aux, aux, s);
    
    //Elegimos un numero aleatorio j tal que tenga los bits suficientes para que
    //p tenga 'bits' bits.
    cout << "Bits aux: " << mpz_sizeinbase(aux, 2) << endl;
    cout << "Bits: " << bits << endl;
    bits_j = bits - mpz_sizeinbase(aux, 2);
    if(bits_j < 1){
    	bits_j = 1;
    }
    mpz_urandomb(j, state, bits_j);
    cout << "j = " << mpz_get_str (nullptr, 10, j) << endl;
    
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
		generate_strong_prime(p, bits/2, state);
		cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
		cout << "tamanio p = " << mpz_sizeinbase(p, 2) << endl;
		cout << "Buscando q..." << endl;
		generate_strong_prime(q, bits/2, state);
		cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
		cout << "Tamanio q = " << mpz_sizeinbase(q, 2) << endl;
		
	}else{
		cout << "Buscando p..." << endl;
		generate_prime(p, bits/2, state);
		cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
		cout << "tamanio p = " << mpz_sizeinbase(p, 2) << endl;
		cout << "Buscando q..." << endl;
		generate_prime(q, bits/2, state);
		cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
		cout << "Tamanio q = " << mpz_sizeinbase(q, 2) << endl;
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
		mpz_urandomb(e, state, mpz_sizeinbase(phi_n, 2)-1);
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

// Factoriza n si p y q son muy cercanos
void ataqueFermat(mpz_t p, mpz_t q, mpz_t n, gmp_randstate_t state){
	mpz_t x, u, v, r;
	mpz_inits(x, u, v, r, nullptr);
	int perf_square = 0;
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		mpz_set(p, n);
		mpz_set_ui(q, 1);
	}
	else{
		//Calculamos el valor inicial sqrt(n)
		//Comprobamos si n es un cuadrado perfecto
		perf_square = mpz_root (x, n, 2);
		if(perf_square != 0){	//Si lo es, damos como solucion la raiz
			mpz_set(p, x);
			mpz_set(q, x);
		}
		
		//u = 2x + 1
		mpz_mul_ui(u, x, 2);
		mpz_add_ui(u, u, 1);
		//v = 1
		mpz_set_ui(v, 1);
		//r = x^2 - n
		mpz_mul(r, x, x);
		mpz_sub(r, r, n);
		
		while(mpz_cmp_si(r,0) != 0){
			if(mpz_cmp_si(r,0) > 0){
				mpz_sub(r,r,v);
				mpz_add_ui(v,v,2);
			}
			else{
				mpz_add(r,r,u);
				mpz_add_ui(u,u,2);
			}
		}
		
		//p = (u+v-2)/2
		mpz_add(p, u, v);
		mpz_sub_ui(p, p, 2);
		mpz_fdiv_q_ui(p, p, 2);
		//q = (u-v)/2
		mpz_sub(q, u, v);
		mpz_fdiv_q_ui(q, q, 2);
	}
	mpz_clears(x, u, v, r, nullptr);
}

// Factoriza n si p y q son muy cercanos
void ataqueKraitchik(mpz_t p, mpz_t q, mpz_t n, gmp_randstate_t state){
	mpz_t x, sq_x;
	mpz_inits(x, sq_x, nullptr);
	int perf_square = 0;
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		mpz_set(p, n);
		mpz_set_ui(q, 1);
	}
	else{
		//Calculamos el valor inicial sqrt(n) (con redondeo hacia arriba)
		//A su vez, comprobamos si n es un cuadrado perfecto
		perf_square = mpz_root (x, n, 2);
		if(perf_square != 0){	//Si lo es, damos como solucion la raiz
			mpz_set(p, x);
			mpz_set(q, x);
		}
		
		mpz_add_ui(x, x, 1);
		while(mpz_cmp(x,n) < 0 && perf_square == 0){	//Mientras x < n, calculamos los cuadrados
			//cout << "x = " << mpz_get_str (nullptr, 10, x) << endl;
			// sq_x = x^2 - n
			mpz_mul(sq_x, x, x);
			mpz_sub(sq_x, sq_x, n);
			//cout << "sq = " << mpz_get_str (nullptr, 10, sq_x) << endl;
			
			perf_square = mpz_root (sq_x, sq_x, 2);	
			//cout << "Perfect square: " << perf_square << endl;
			if(perf_square != 0){	// Si sq_x es un cuadrado perfecto
				mpz_add(p,x,sq_x);	// p = x+y
				mpz_sub(q,x,sq_x);	// q = x-y
			}
			mpz_add_ui(x,x,1);
		}
	}
	mpz_clears(x, sq_x, nullptr);
}

//Funcion con comportamiento pseudoaleatorio que emplea el metodo rho de Pollard
void paso_rhoPollard(mpz_t x, mpz_t n){
	//Funcion: f(x) = (x^2 + 1) mod n
	mpz_mul(x, x, x);
	mpz_add_ui(x, x, 1);
	mpz_fdiv_r(x, x, n);
}

// Factoriza n si p y q son muy cercanos
bool rhoPollard(mpz_t p, mpz_t q, mpz_t n, gmp_randstate_t state){
	mpz_t t, l, dif, mcd;
	mpz_inits(t, l, dif, mcd, nullptr);
	unsigned long int init = 2;
	bool exito = true;
	
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		mpz_set(p, n);
		mpz_set_ui(q, 1);
	}
	else{
		mpz_set_ui(t, init);
		mpz_set_ui(l, init);
		mpz_set_ui(mcd, 1);
		
		while(mpz_cmp_ui(mcd, 1) == 0){
			//Avanzamos la tortuga un paso
			paso_rhoPollard(t, n);
			//Avanzamos la liebre dos pasos
			paso_rhoPollard(l, n);
			paso_rhoPollard(l, n);
			//Calculamos el mcd de la diferencia t - l con n
			mpz_sub(dif, t, l);
			if(mpz_cmp_ui(dif, 0) < 0){
				mpz_mul_si(dif, dif, -1);
			}
			cout << "t = " << mpz_get_str (nullptr, 10, t) << endl;
			cout << "l = " << mpz_get_str (nullptr, 10, l) << endl;
			cout << "|t - l| = " << mpz_get_str (nullptr, 10, dif) << endl << flush;
			mpz_gcd (mcd, n, dif);
		}
		mpz_set(p, mcd);
		mpz_fdiv_q(q, n, mcd);	
		
		if(mpz_cmp_ui(q, 1) == 0){
			exito = false;
		}
	}
	mpz_clears(t, l, dif, mcd, nullptr);
	
	return exito;
}


struct Punto{
	mpz_t x;
	mpz_t y;
};

//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
int sumaCurvaEliptica(Punto &suma, const Punto s1, const Punto s2, const mpz_t a, 
					  const mpz_t b, const mpz_t p, mpz_t inv){
	mpz_t aux, lambda;
	Punto s3;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente
	mpz_inits(aux, lambda, s3.x, s3.y, nullptr);
	//cout << "Suma curva eliptica" << endl;
	
	//cout << "s1 = [" << mpz_get_str (nullptr, 10, s1.x) << ", " 
	//				  << mpz_get_str (nullptr, 10, s1.y) << "]" << endl;		  
	//cout << "s2 = [" << mpz_get_str (nullptr, 10, s2.x) << ", " 
	//				  << mpz_get_str (nullptr, 10, s2.y) << "]" << endl;
	//cout << "p = " << mpz_get_str (nullptr, 10, s2.x) << endl;
	//Si alguno de los puntos es O, devolvemos el otro punto como solucion
	if(mpz_cmp_ui(s1.x, 0) < 0){
		mpz_set(s3.x, s2.x);
		mpz_set(s3.y, s2.y);
		//cout << "Sumando 1 es \"O\": devolviendo s2" << endl;
	}
	else{
		if(mpz_cmp_ui(s2.x, 0) < 0){
			mpz_set(s3.x, s1.x);
			mpz_set(s3.y, s1.y);
			//cout << "Sumando 2 es \"O\": devolviendo s1" << endl;
		}
		else{
			//Si tenemos dos puntos de la forma (x,y) y (x,-y), devolvemos O
			
			mpz_add(aux, s1.y, s2.y);
			//cout << "Suma curva eliptica sin O" << endl;
			mpz_mod(aux, aux, p);
			//cout << "Suma curva eliptica sin O" << endl;
			if(mpz_cmp(s1.x, s2.x) == 0 && mpz_cmp_ui(aux, 0) == 0){
				mpz_set_si(s3.x, -1);
				mpz_set_si(s3.y, -1);
				//cout << "Puntos opuestos: devolviendo \"O\"" << endl;
			}
			//En otro caso
			else{
				//Calculamos lambda
				//Si ambos puntos son iguales
				if(mpz_cmp(s1.x, s2.x) == 0 && mpz_cmp(s1.y, s2.y) == 0){
					//cout << "Puntos iguales: calculando lambda" << endl;
					//lambda = (3x^2 + a)/2y
					mpz_mul(lambda, s1.x, s1.x);
					mpz_mul_ui(lambda, lambda, 3);
					mpz_add(lambda, lambda, a);
					mpz_mul_ui(inv, s1.y, 2);
					exito = mpz_invert(aux, inv, p);	//Si no hay inversa, devuelve 0
					mpz_mul(lambda, lambda, aux);
					mpz_mod(lambda, lambda, p);
				}
				//Si ambos numeros son diferentes
				else{
					//cout << "Puntos diferentes: calculando lambda" << endl;
					//lambda = (x1 - x2)/(y1 - y2)
					mpz_sub(lambda, s1.y, s2.y);
					mpz_sub(inv, s1.x, s2.x);
					exito = mpz_invert(aux, inv, p);	//Si no hay inversa, devuelve 0
					//if(exito == 0)
						//cout << "No existe el inverso de " << inv << " en Z_" << p << endl;
					
					mpz_mul(lambda, lambda, aux);
					mpz_mod(lambda, lambda, p);
				}
				
				//Calculamos la coordenada x de la suma
				//x3 = lambda^2 - x1 - x2 (mod p)
				mpz_mul(s3.x, lambda, lambda);
				mpz_sub(s3.x, s3.x, s1.x);
				mpz_sub(s3.x, s3.x, s2.x);
				mpz_mod(s3.x, s3.x, p);
				
				//Calculamos la coordenada y
				//y3 = lambda*(x_1-x_3) - y1 (mod p)
				mpz_sub(s3.y, s1.x, s3.x);
				mpz_mul(s3.y, s3.y, lambda);
				mpz_sub(s3.y, s3.y, s1.y);
				mpz_mod(s3.y, s3.y, p);
			}
		}
	}
	//Asignamos el resultado a la salida
	mpz_set(suma.x, s3.x);
	mpz_set(suma.y, s3.y);
	
	mpz_clears(aux, lambda, s3.x, s3.y, nullptr);
	return exito;
}

int multiplicacionCurvaEliptica(Punto &mul, const Punto f1, const mpz_t f2, const mpz_t a, 
								const mpz_t b, const mpz_t p, mpz_t inv){
	Punto p1, p2;
	mpz_t aux, lambda, f2_2;
	mp_limb_t f2_actual;
	int size_f2, bit_actual, resto;
	int size_limb = sizeof(mp_limb_t)*8;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente en las sumas
	int i = 0;
	bool primer_uno_encontrado = false;
	mpz_inits(aux, lambda, p1.x, p1.y, p2.x, p2.y, f2_2, nullptr);
	//cout << "Inicio multiplicacion" << endl;
	
	//Inicializamos la multiplicacion en "O"
	mpz_set_si(mul.x, -1);
	mpz_set_si(mul.y, -1);
	//Guardamos el punto a multiplicar P en p2
	mpz_set(p2.x, f1.x);
	mpz_set(p2.y, f1.y);
	mpz_set(f2_2, f2);
	//cout << "p2.x = " << mpz_get_str (nullptr, 10, p2.x) << endl;
	//cout << "p2.y = " << mpz_get_str (nullptr, 10, p2.y) << endl;
	
	while(mpz_cmp_ui(f2_2,0) > 0 && exito != 0){
		//cout << "Iteracion " << i << endl;
		//Calculamos el siguiente bit
		//cout << "f2_2 = " << mpz_get_str (nullptr, 2, f2_2) << endl;
		resto = mpz_fdiv_qr_ui(f2_2, aux, f2_2, 2);
		//cout << "Resto: " << resto << endl;
		if(mpz_cmp_ui(aux, 1) == 0){
			exito = sumaCurvaEliptica(mul, mul, p2, a, b, p, inv);
			//cout << "Exito 1: " << exito << endl;
		}
		//Calculamos la siguiente potencia de 2 de P
		if(exito != 0){
			exito = sumaCurvaEliptica(p1, p2, p2, a, b, p, inv);
			mpz_set(p2.x, p1.x);
			mpz_set(p2.y, p1.y);
			//cout << "Exito 2: " << exito << endl;
		}
		i++;
	}
	
	mpz_clears(aux, lambda, p1.x, p1.y, p2.x, p2.y, nullptr);

	return exito;
}

//Calcula el numero mas grande tal que sea k-potencia-suave
void maxKPotenciaSuave(mpz_t kps, const mpz_t k){
	mpz_t p_actual, p_potencia;
	mpz_inits (p_actual, p_potencia, nullptr);
	
	mpz_set_ui(kps, 1);
	mpz_set_ui(p_actual, 2);
	
	//Mientras el primo actual no supere el valor de K
	while(mpz_cmp(p_actual, k) < 0){
		mpz_set(p_potencia, p_actual);
		//Calculamos el valor de p_actual^e tal que no supere K
		while(mpz_cmp(p_potencia, k) < 0){
			mpz_mul(p_potencia,p_potencia,p_actual);	//p_potencia = p_actual^{n+1}
		}
		mpz_divexact(p_potencia, p_potencia, p_actual);	//Dividimos para que p_potencia < k
		//Multiplicamos dicho valor en el acumulador
		mpz_mul(kps, kps, p_potencia);
		//Pasamos al siguiente primo
		mpz_nextprime(p_actual, p_actual);
	}
	mpz_clears(p_actual, p_potencia, nullptr);
}

bool factorizacionCurvasElipticas(mpz_t p, mpz_t q, const mpz_t n, gmp_randstate_t state, const mpz_t k, const unsigned int att){
	Punto Q;
	mpz_t a, b, aux, inv, f2;
	mpz_inits(a, b, Q.x, Q.y, aux, inv, f2, nullptr);
	mp_limb_t random;
	int encontrado_inv = 1;	//Comprueba si el calculo de la inversa se hizo correctamente en las sumas
	bool mul_O = false; 	//Indica si hemos llegado a obtener O como resultado de la multiplicacion
	bool exito = false;
	
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		mpz_set(p, n);
		mpz_set_ui(q, 1);
		exito = true;
	}
	else{
		for(int i = 0; i < att && !exito; i++){
			//Escogemos una pseudo curva eliptica eligiendo a, b y definiendola en Z_n
			mpz_urandomm (a, state, n);
			mpz_urandomm (b, state, n);
			//Escogemos un punto Q aleatorio de la curva
			mpz_urandomm (Q.x, state, n);	//Elegimos x aleatorio
			//Calculamos y como x^3 + ax + b (mod n)
			mpz_mul(aux, Q.x, Q.x);
			mpz_mul(Q.y, aux, Q.x);
			mpz_addmul(Q.y, a, Q.x);
			mpz_add(Q.y, Q.y, b);
			mpz_mod(Q.y, Q.y, n);
			
			mpz_set_ui(f2, 2);
			
			//Usamos el factorial
			/*
			while(encontrado_inv != 0 && !mul_O && mpz_cmp(f2, k) < 0){
				encontrado_inv = multiplicacionCurvaEliptica(Q, Q, f2, a, b, n, inv);
				mpz_add_ui(f2, f2, 1);
				if(mpz_cmp_ui(Q.x, 0) < 0){
					mul_O = true;
				}
			}
			*/
			//Calculamos L = maximo numero K-potencia-uniforme. 
			maxKPotenciaSuave(f2, k);
			//cout << "hola 1 " << endl << flush;
			//Realizamos la multiplicacion
			encontrado_inv = multiplicacionCurvaEliptica(Q, Q, f2, a, b, n, inv);
			//cout << "hola 2 " << endl<< flush;
			//Si hemos encontrado un elemento no invertible, obtenemos un factor de n
			if(encontrado_inv == 0){
				//Calculamos gcd
				mpz_gcd (p, inv, n);
				mpz_fdiv_q(q, n, p);
				exito = true;
			}
		}
	}
	mpz_clears(a, b, Q.x, Q.y, aux, inv, f2, nullptr);
	return exito;
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
	bool exito;
	int att = 10;
	mpz_t n, p, q, k;
	mpz_inits(p, q, n, k, nullptr);
	mpz_set_ui(n, 265905204186593);
	mpz_set_ui(k, 50);
	exito = factorizacionCurvasElipticas(p, q, n, state, k, att);
	if(exito){
		cout << "n = " << mpz_get_str (nullptr, 10, n) << endl;
		cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
		cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
	}
	else{
		cout << "No se ha conseguido factorizar el numero :(" << endl;
	}
	mpz_clears(n, p, q, k, nullptr);
	
	
    return 0;
}

