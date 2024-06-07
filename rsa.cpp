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
bool millerRabin(mpz_t p, int it, gmp_randstate_t state){
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
void ataqueFermat(mpz_t n, mpz_t p, mpz_t q){
	mpz_t x, u, v, r;
	mpz_inits(x, u, v, r, nullptr);
	int perf_square = 0;
	
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
	
	mpz_clears(x, u, v, r, nullptr);
}

// Factoriza n si p y q son muy cercanos
void ataqueKraitchik(mpz_t n, mpz_t p, mpz_t q){
	mpz_t x, sq_x;
	mpz_inits(x, sq_x, nullptr);
	int perf_square = 0;
	
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
bool rhoPollard(mpz_t n, mpz_t p, mpz_t q){
	mpz_t t, l, dif, mcd;
	mpz_inits(t, l, dif, mcd, nullptr);
	unsigned long int init = 2;
	bool exito = true;
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
	mpz_clears(t, l, dif, mcd, nullptr);
	
	return exito;
}


struct Punto{
	mpz_t x;
	mpz_t y;
};

//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
void sumaCurvaEliptica(Punto suma, Punto s1, Punto s2, mpz_t a, mpz_t b, mpz_t p){
	mpz_t aux, lambda;
	Punto s3;
	mpz_inits(aux, lambda, s3.x, s3.y, nullptr);
	
	//Si alguno de los puntos es O, devolvemos el otro punto como solucion
	if(mpz_cmp_ui(s1.x, 0) < 0){
		mpz_set(suma.x, s2.x);
		mpz_set(suma.y, s2.y);
	}
	else{
		if(mpz_cmp_ui(s2.x, 0) < 0){
			mpz_set(suma.x, s1.x);
			mpz_set(suma.y, s1.y);
		}
		else{
			//Si tenemos dos puntos de la forma (x,y) y (x,-y), devolvemos O
			mpz_add(aux, s1.y, s2.y);
			if(mpz_cmp(s1.x, s2.x) == 0 && mpz_cmp_ui(aux, 0) == 0){
				mpz_set_si(suma.x, -1);
				mpz_set_si(suma.y, -1);
			}
			//En otro caso
			else{
				//Calculamos lambda
				//Si ambos puntos son iguales
				if(mpz_cmp(s1.x, s2.x) == 0 && mpz_cmp(s1.y, s2.y) == 0){
					//lambda = (3x^2 + a)/2y
					mpz_mul(lambda, s1.x, s1.x);
					mpz_mul_ui(lambda, lambda, 3);
					mpz_add(lambda, lambda, a);
					mpz_mul_ui(aux, s1.y, 2);
					mpz_invert(aux, aux, p);
					mpz_mul(lambda, lambda, aux);
				}
				//Si ambos numeros son diferentes
				else{
					//lambda = (x1 - x2)/(y1 - y2)
					mpz_sub(lambda, s1.y, s2.y);
					mpz_sub(aux, s1.x, s2.x);
					mpz_invert(aux, aux, p);
					mpz_mul(lambda, lambda, aux);
				}
				
				//Calculamos la coordenada x de la suma
				//x3 = lambda^2 - x1 - x2 (mod p)
				mpz_mul(s3.x, lambda, lambda);
				mpz_sub(s3.x, s3.x, s1.x);
				mpz_sub(s3.x, s3.x, s2.x);
				mpz_fdiv_r(suma.x, suma.x, p);
				
				//Calculamos la coordenada y
				//y3 = lambda*(x_1-x_3) - y1 (mod p)
				mpz_sub(s3.y, s1.x, s3.x);
				mpz_mul(s3.y, s3.y, lambda);
				mpz_sub(s3.y, s3.y, s1.y);
			}
		}
	}
	//Asignamos el resultado a la salida
	mpz_set(suma.x, s3.x);
	mpz_set(suma.y, s3.y);
	
	mpz_clears(aux, lambda, s3.x, s3.y, nullptr);
	
}

void multiplicacionCurvaEliptica(Punto mul, Punto f1, mpz_t f2, mpz_t a, mpz_t b, mpz_t p){
	Punto p1, p2;
	mpz_t aux, lambda;
	mp_limb_t f2_actual;
	int size_f2, bit_actual;
	int size_limb = sizeof(mp_limb_t)*8;
	bool primer_uno_encontrado = false;
	mpz_inits(aux, lambda, p1.x, p1.y, p2.x, p2.y, nullptr);
	
	//Guardamos el punto inicial f1 en p1 y calculamos p2=2*p1=p1+p1
	mpz_set(mul.x, f1.x);
	mpz_set(mul.y, f1.y);
	sumaCurvaEliptica(p2, p1, p1, a, b, p);
	
	//Obtenemos la representacion en binario de f2
	size_f2 = mpz_size(f2);
	for(int i = size_f2 - 1; i >=0 ; i--){
		f2_actual = mpz_getlimbn(f2, i);
		//Iteramos sobre cada bit de f2_actual
		for(int j = 0; j < size_limb; j--){
			//Nos quedamos con el bit size_limb - j - 1
			bit_actual = (bit_actual << j) >> (size_limb-1);
			//Si no hemos encontrado aun el primer bit a 1, no hacemos nada
			if(!primer_uno_encontrado){
				if(bit_actual > 0){
					primer_uno_encontrado = true;
				}
			}
			else{
				//Si el bit actual es 1
				if(bit_actual > 0){
					sumaCurvaEliptica(mul, mul, p2, a, b, p);	//p1 = p1 + p2
					sumaCurvaEliptica(p2, p2, p2, a, b, p);		//p2 = 2p2
				}
				//Si el bit actual es 0
				else{
					sumaCurvaEliptica(p2, mul, p2, a, b, p);	//p2 = p1 + p2
					sumaCurvaEliptica(mul, mul, mul, a, b, p);	//p1 = 2p1
				}
			}
		}
	}
	mpz_clears(aux, lambda, p1.x, p1.y, p2.x, p2.y, nullptr);
}

void factorizacionCurvasElipticas(mpz_t n, mpz_t p, mpz_t q){
	//Escogemos una pseudo curva eliptica eligiendo a, b y definiendola en Z_p
	//Escogemos un punto Q aleatorio de la curva
	//Calculamos L = maximo numero K-potencia-uniforme. TODO Funcion que devuelva dicho valor para cierto K
	//Calculamos L*Q
	//Si mientras se calcula se da que l*Q = O --> l | n --> l | p y l | q (bajar la cota)
	//Si llegamos a l*Q = A+B: l | p --> En curva eliptica mod p, obtenemos O --> lambda es infinito -->
	//x_1 - x_2 = 0 (mod p) --> x_1 - x_2 no tiene inverso, por lo que no se puede calcular lambda -->
	//mcd(x_1 - x_2, N) = p.
	//Si no llegamos a ninguno de los anteriores --> Subir cota K / Cambiar curva.
	
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
	// Prueba Fermat
	mpz_t n, p, q;
	mpz_inits(n, p, q, nullptr);
	mpz_set_ui(n, 63205667035);
	rhoPollard(n, p, q);
	cout << "n = " << mpz_get_str (nullptr, 10, n) << endl;
	cout << "p = " << mpz_get_str (nullptr, 10, p) << endl;
	cout << "q = " << mpz_get_str (nullptr, 10, q) << endl;
	mpz_clears(n, p, q, nullptr);
	
	
    return 0;
}

