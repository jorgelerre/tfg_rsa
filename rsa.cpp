#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

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
bool strongPseudoprime(const mpz_class p, const mpz_class base, bool debug = false){
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
	unsigned long int a;
	mpz_class b;
	
	//Realizamos el test de las divisiones sucesivas con los primeros primos (<8 bits)
	const int NUM_PRIMES = 54;
	const unsigned int primes[NUM_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
			   59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
	 		   157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};
	
	if(p > 252){
		for(int i = 0; i < NUM_PRIMES && no_fail; i++){
			b = p % primes[i];
			if(b == 0)
				no_fail = false;
		}
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
mpz_class generate_prime(int bits, gmp_randstate_t state) {
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
mpz_class generate_strong_prime(int bits, gmp_randstate_t state) {
    mpz_class p,p_0,r,s,t,i,j,aux;
    int bits_j, bits_1, bits_2;
    //cout << "Buscando s,t..." << endl;
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
		//cout << "s = " << mpz_get_str (nullptr, 10, s) << endl;
		//cout << "t = " << mpz_get_str (nullptr, 10, t) << endl;
		
		//Elegimos un numero aleatorio i
		mpz_urandomb(i.get_mpz_t(), state, log2(bits));
		//cout << "i = " << mpz_get_str (nullptr, 10, i) << endl;
		//Calculamos r: r sea primo
		do{
			i++;
			r = 2*i*t + 1;
		}while(!millerRabin(r,10,state));
		//cout << "r = " << mpz_get_str (nullptr, 10, r) << endl;
		//cout << "tamanio r = " << mpz_sizeinbase(r, 2) << endl;
		//Calculamos p_0
		//p_0 = ((2*s^{r-2}) % r) *s - 1
		p_0 = r - 2;
		mpz_powm_sec(p_0.get_mpz_t(), s.get_mpz_t(), p_0.get_mpz_t(), r.get_mpz_t());	//p_0 = s^(p_0) (mod r)
		//cout << "p_02 = " << mpz_get_str (nullptr, 10, p_0) << endl;
		p_0 = ((2*p_0) % r)*s - 1;

		//if(!mpz_odd_p(p_0))
			//cout << "p_0 no es impar: repetimos el proceso" << endl;
	}while(p_0 % 2 == 0);	//p_0 debe ser impar, si no, repetimos el proceso
	/*
	cout << "tamanio s = " << mpz_sizeinbase(s, 2) << endl;
	cout << "tamanio t = " << mpz_sizeinbase(t, 2) << endl;
	cout << "tamanio p_0 = " << mpz_sizeinbase(p_0, 2) << endl;
    */
    //aux = 2rs
    aux = 2*r*s;
    
    //Elegimos un numero aleatorio j tal que tenga los bits suficientes para que
    //p tenga 'bits' bits.
    //cout << "Bits aux: " << mpz_sizeinbase(aux.get_mpz_t(), 2) << endl;
    //cout << "Bits: " << bits << endl;
    bits_j = bits - mpz_sizeinbase(aux.get_mpz_t(), 2);
    if(bits_j < 1){
    	bits_j = 1;
    }
    mpz_urandomb(j.get_mpz_t(), state, bits_j);
    //cout << "j = " << mpz_get_str (nullptr, 10, j) << endl;
    
    //p = 2rsj + p_0
    p = aux*j + p_0; 
    /*
    cout << "aux = " << aux << endl;
    cout << "j = " << j << endl;
    cout << "p = " << p << endl;
    */
    //Mientras p no sea primo, le sumamos aux
    while(!millerRabin(p,10,state)){
    	p += aux;
    	//cout << "p = " << p << endl;
    }
    
    return p;
}

// Función principal para generar claves RSA
void generate_rsa_key(int bits, mpz_class &n, mpz_class &e, mpz_class &d, gmp_randstate_t state, bool strong_prime = true) {
	mpz_class p, q, phi_n, mcd;
	
	unsigned int e_default = 65537;
	//Elegimos dos primos p y q aleatoriamente
	if(strong_prime){
		cout << "Buscando p..." << endl;
		p = generate_strong_prime(bits/2, state);
		cout << "p = " << p << endl;
		cout << "tamanio p = " << mpz_sizeinbase(p.get_mpz_t(), 2) << endl;
		cout << "Buscando q..." << endl;
		q = generate_strong_prime(bits/2, state);
		cout << "q = " << q << endl;
		cout << "Tamanio q = " << mpz_sizeinbase(q.get_mpz_t(), 2) << endl;
	}else{
		cout << "Buscando p..." << endl;
		p = generate_prime(bits/2, state);
		cout << "p = " << p << endl;
		cout << "tamanio p = " << mpz_sizeinbase(p.get_mpz_t(), 2) << endl;
		cout << "Buscando q..." << endl;
		q = generate_prime(bits/2, state);
		cout << "q = " << q << endl;
		cout << "Tamanio q = " << mpz_sizeinbase(q.get_mpz_t(), 2) << endl;
	}
	//Calculamos n y phi(n)
	cout << "Calculando n y phi(n)..." << endl;
	n = p*q;			
	phi_n = (p-1)*(q-1);

	cout << "n = " << n << endl;
	cout << "phi(n) = " << phi_n << endl;

	cout << "Calculando e..." << endl;
	//Elegimos el exponente de cifrado e = 2^16-1 = 65.535 si es < phi_n
	if(phi_n > e_default){
		cout << "Asignando e_default..." << endl;
		e = e_default;
	}
	//En otro caso, lo elegimos aleatoriamente (un bit mas pequeño que phi_n)
	else{
		cout << "Asignando e aleatorio..." << endl;
		mpz_urandomb(e.get_mpz_t(), state, mpz_sizeinbase(phi_n.get_mpz_t(), 10)-1);
		if(e % 2 == 0){
			e--;
		}
	}
	
	mpz_gcd(mcd.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
	cout << "mcd(e,phi_n) = " << mcd << endl;
	//Mientras e y phi_n no sean coprimos
	while(mcd > 1 || e < 3 || e > phi_n){
		if(mcd > 1)
			cout << "Asignando e aleatorio (eran coprimos)..." << endl;
		if(e < 3)
			cout << "Asignando e aleatorio (e era menor que 3)..." << endl;
		if(e > phi_n){
			cout << "Asignando e aleatorio (e era mayor que phi(n))..." << endl;
			cout << "phi(n) = " << phi_n << endl;
		}
		mpz_urandomb(e.get_mpz_t(), state, mpz_sizeinbase(phi_n.get_mpz_t(), 2)-1);
		if(e % 2 == 0){
			e--;
		}
		
		mpz_gcd(mcd.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
		cout << "e      = " << e << endl;
		cout << "mcd(e,phi_n) = " << mcd << endl;
	}
	
	cout << "Calculando d..." << endl;
	
	//Calculamos el inverso de e en Z_phi(n) = d
	mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi_n.get_mpz_t());
	
	//TODO: Comprobar ataque Wiener
	
}

mpz_class cifra_RSA(mpz_class m, mpz_class e, mpz_class n){
	mpz_class m_cifrado;
	//Hacemos la operacion m_cifrado = m^e (mod n)
	mpz_powm(m_cifrado.get_mpz_t(), m.get_mpz_t(), e.get_mpz_t(), n.get_mpz_t());
	return m_cifrado;
}

mpz_class descifra_RSA(mpz_class m, mpz_class e, mpz_class n){
	return cifra_RSA(m, e, n);
}

// Factoriza n si p y q son muy cercanos
mpz_class ataqueFermat(mpz_class n, gmp_randstate_t state){
	mpz_class p, x, u, v, r;
	
	int perf_square = 0;
	//Comprobamos que n no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		//Calculamos el valor inicial sqrt(n)
		//Comprobamos si n es un cuadrado perfecto
		perf_square = mpz_root(x.get_mpz_t(), n.get_mpz_t(), 2);
		if(perf_square != 0){	//Si lo es, damos como solucion la raiz
			p = x;
		}
		//Si no, ejecutamos el ataque de Fermat
		else{
			u = 2*x + 1;
			v = 1;
			//r = x^2 - n
			r = x*x - n;
			
			while(r != 0){
				if(r > 0){
					r = r - v;
					v = v + 2;
				}
				else{
					r = r + u;
					u = u + 2;
				}
			}
			
			//p = (u+v-2)/2
			p = (u+v-2)/2;
		}
	}
	
	return p;
}


mpz_class ataqueKraitchik(mpz_class n, gmp_randstate_t state){
	mpz_class x, sq_x, p;
	
	int perf_square = 0;
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		//Calculamos el valor inicial sqrt(n) (con redondeo hacia arriba)
		//A su vez, comprobamos si n es un cuadrado perfecto
		perf_square = mpz_root(x.get_mpz_t(), n.get_mpz_t(), 2);
		if(perf_square != 0){	//Si lo es, damos como solucion la raiz
			p = x;
		}
		
		while(x < n && perf_square == 0){	//Mientras x < n y no haya solucion, calculamos los cuadrados
			x++;
			//cout << "x = " << x << endl;
			sq_x = x*x - n;
			//cout << "sq = " << sq_x << endl;
			
			perf_square = mpz_root(sq_x.get_mpz_t(), sq_x.get_mpz_t(), 2);	
			//cout << "Perfect square: " << perf_square << endl;
			if(perf_square != 0){	// Si sq_x es un cuadrado perfecto
				p = x + sq_x;
			}
		}
	}
	return p;
}

//Funcion con comportamiento pseudoaleatorio que emplea el metodo rho de Pollard
mpz_class paso_rhoPollard(mpz_class x, mpz_class n){
	//Funcion: f(x) = (x^2 + 1) mod n
	return (x*x + 1) % n;
}

// Factoriza n si p y q son muy cercanos
mpz_class rhoPollard(mpz_class n, gmp_randstate_t state){
	mpz_class p, t, l, dif, mcd = 1;
	
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		//Tortuga y liebre empiezan en el mismo estado (2)
		t = 2;
		l = 2;
		
		while(mcd == 1){
			//Avanzamos la tortuga un paso
			t = paso_rhoPollard(t, n);
			//Avanzamos la liebre dos pasos
			l = paso_rhoPollard(paso_rhoPollard(l, n), n);

			//Calculamos el mcd de la diferencia t - l con n
			dif = t - l;
			if(dif < 0){
				dif *= -1;
			}
			cout << "t = " << t << endl;
			cout << "l = " << l << endl;
			cout << "|t - l| = " << dif << endl << flush;
			mpz_gcd (mcd.get_mpz_t(), n.get_mpz_t(), dif.get_mpz_t());
		}
		p = mcd;
		
	}
	
	return p;
}


struct Punto{
	mpz_class x;
	mpz_class y;
};

//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
Punto sumaCurvaEliptica(const Punto s1, const Punto s2, const mpz_class a, 
					  const mpz_class b, const mpz_class p, mpz_class &inv){
	Punto s3;
	mpz_class aux, lambda;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente
	
	//cout << "Suma curva eliptica" << endl;
	//cout << "s1 = [" <<  s1.x << ", " 
	//				  <<  s1.y << "]" << endl;		  
	//cout << "s2 = [" << s2.x << ", " 
	//				  << s2.y << "]" << endl;
	//cout << "p = " << p << endl;
	//Si alguno de los puntos es O, devolvemos el otro punto como solucion
	if(s1.x == 0){
		s3.x = s2.x;
		s3.y = s2.y;
		//cout << "Sumando 1 es \"O\": devolviendo s2" << endl;
	}
	else{
		if(s2.x < 0){
			s3.x = s1.x;
			s3.y = s1.y;
			//cout << "Sumando 2 es \"O\": devolviendo s1" << endl;
		}
		else{
			//Si tenemos dos puntos de la forma (x,y) y (x,-y), devolvemos O
			aux = s1.y + s2.y;
			//cout << "Suma curva eliptica sin O" << endl;
			mpz_mod(aux.get_mpz_t(), aux.get_mpz_t(), p.get_mpz_t());
			//cout << "Suma curva eliptica sin O" << endl;
			if(s1.x == s2.x && aux == 0){
				s3.x = -1;
				s3.y = -1;
				//cout << "Puntos opuestos: devolviendo \"O\"" << endl;
			}
			//En otro caso
			else{
				//Calculamos lambda
				//Si ambos puntos son iguales
				if(s1.x == s2.x && s1.y == s2.y){
					//cout << "Puntos iguales: calculando lambda" << endl;
					//lambda = (3x^2 + a)/2y
					lambda = 3*s1.x*s1.x + a;
					inv = 2*s1.y;
					//Si no hay inversa, devuelve 0
					exito = mpz_invert(aux.get_mpz_t(), inv.get_mpz_t(), p.get_mpz_t());	
					if(exito != 0){
						lambda = lambda * aux;
						mpz_mod(lambda.get_mpz_t(), lambda.get_mpz_t(), p.get_mpz_t());
					}
				}
				//Si ambos numeros son diferentes
				else{
					//cout << "Puntos diferentes: calculando lambda" << endl;
					//lambda = (y1 - y2)/(x1 - x2)
					lambda = s1.y - s2.y;
					inv = s1.x - s2.x;
					//Si no hay inversa, devuelve 0
					exito = mpz_invert(aux.get_mpz_t(), inv.get_mpz_t(), p.get_mpz_t());	
					if(exito != 0){
						lambda = lambda*aux;	
						mpz_mod(lambda.get_mpz_t(), lambda.get_mpz_t(), p.get_mpz_t());
					}
				}
				if(exito != 0){
					//Calculamos la coordenada x de la suma
					//x3 = lambda^2 - x1 - x2 (mod p)
					s3.x = lambda*lambda - s1.x - s2.x;
					mpz_mod(s3.x.get_mpz_t(), s3.x.get_mpz_t(), p.get_mpz_t());
					
					//Calculamos la coordenada y
					//y3 = lambda*(x_1-x_3) - y1 (mod p)
					s3.y = lambda*(s1.x - s3.x) - s1.y;
					mpz_mod(s3.y.get_mpz_t(), s3.y.get_mpz_t(), p.get_mpz_t());
				}
				else{
					//Punto no definido
					s3.x = -2;
					s3.y = -2;
				}
			}
		}
	}
	
	return s3;
}

Punto multiplicacionCurvaEliptica(const Punto f1, const mpz_class f2, const mpz_class a, 
								const mpz_class b, const mpz_class p, mpz_class &inv){
	Punto p2, mul;
	mpz_class r, lambda, f2_2;
	mp_limb_t f2_actual;
	int size_f2, bit_actual, resto;
	int size_limb = sizeof(mp_limb_t)*8;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente en las sumas
	int i = 0;
	bool primer_uno_encontrado = false;
	
	//cout << "Inicio multiplicacion" << endl;
	
	//Inicializamos la multiplicacion en "O"
	mul.x = -1;
	mul.y = -1;
	
	//Guardamos el punto a multiplicar P en p2
	p2.x = f1.x;
	p2.y = f1.y;
	f2_2 = f2;
	//cout << "p2.x = " << p2.x << endl;
	//cout << "p2.y = " << p2.y << endl;
	
	while(f2_2 > 0 && exito != 0){
		//cout << "Iteracion " << i << endl;
		//Calculamos el siguiente bit
		//cout << "f2_2 = " << f2_2 << endl;
		resto = mpz_fdiv_qr_ui(f2_2.get_mpz_t(), r.get_mpz_t(), f2_2.get_mpz_t(), 2);
		//cout << "Resto: " << resto << endl;
		if(r == 1){
			mul = sumaCurvaEliptica(mul, p2, a, b, p, inv);
			if(mul.x == -2){
				exito = 0;
			}
			//cout << "Exito 1: " << exito << endl;
		}
		//Calculamos la siguiente potencia de 2 de P
		if(exito != 0){
			p2 = sumaCurvaEliptica(p2, p2, a, b, p, inv);
			if(p2.x == -2){
				exito = 0;
			}
			//cout << "Exito 2: " << exito << endl;
		}
		i++;
	}
	//Devolvemos un punto no existente si falla alguna suma
	if(exito == 0){
		mul.x = -2;
		mul.y = -2;
	}

	return mul;
}

//Calcula el numero mas grande tal que sea k-potencia-suave
mpz_class maxKPotenciaSuave(const mpz_class k){
	mpz_class kps, p_actual, p_potencia;
	
	kps = 1;
	p_actual = 2;
	
	//Mientras el primo actual no supere el valor de K
	while(p_actual < k){
		p_potencia = p_actual;
		//Calculamos el valor de p_actual^e tal que no supere K
		while(p_potencia < k){
			p_potencia = p_potencia * p_actual;	//p_potencia = p_actual^{n+1}
		}
		//Dividimos para que p_potencia < k
		mpz_divexact(p_potencia.get_mpz_t(), p_potencia.get_mpz_t(), p_actual.get_mpz_t());	
		//Multiplicamos dicho valor en el acumulador
		kps = kps * p_potencia;
		//Pasamos al siguiente primo
		mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
	}
	
	return kps;
}

mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k, const unsigned int att){
	Punto Q, M;
	mpz_class p, a, b, aux, inv, L;
	
	mp_limb_t random;
	int encontrado_inv = 1;	//Comprueba si el calculo de la inversa se hizo correctamente en las sumas
	bool mul_O = false; 	//Indica si hemos llegado a obtener O como resultado de la multiplicacion
	bool exito = false;
	
	//Comprobamos que n no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		for(int i = 0; i < att && !exito; i++){
			//Escogemos una pseudo curva eliptica eligiendo a, b y definiendola en Z_n
			mpz_urandomm (a.get_mpz_t(), state, n.get_mpz_t());
			mpz_urandomm (b.get_mpz_t(), state, n.get_mpz_t());
			//Escogemos un punto Q aleatorio de la curva
			mpz_urandomm (Q.x.get_mpz_t(), state, n.get_mpz_t());	//Elegimos x aleatorio
			//Calculamos y como x^3 + ax + b (mod n)
			Q.y = Q.x * Q.x * Q.x + a*Q.x + b;
			mpz_mod(Q.y.get_mpz_t(), Q.y.get_mpz_t(), n.get_mpz_t());
			
			//Calculamos L = maximo numero K-potencia-uniforme. 
			L = maxKPotenciaSuave(k);
			//cout << "hola 1 " << endl << flush;
			//Realizamos la multiplicacion
			M = multiplicacionCurvaEliptica(Q, L, a, b, n, inv);
			//cout << "hola 2 " << endl<< flush;
			//Si hemos encontrado un elemento no invertible, obtenemos un factor de n
			if(M.x == -2){
				//Calculamos gcd
				mpz_gcd (p.get_mpz_t(), inv.get_mpz_t(), n.get_mpz_t());
				exito = true;
			}
		}
	}
	//Si no hemos conseguido factorizar, devolvemos p = n
	if(!exito){
		p = n;
	}
	
	return p;
}


// Función para encontrar una raíz cuadrada de a módulo p usando el algoritmo de Tonelli-Shanks
void sqrt_mod(mpz_class &sqrt_a, const mpz_class &a, const mpz_class &n){
	mpz_class a2, b, q = n-1, s = 0, z, aux;
	mpz_class M, c, t, R;
	//cout << "a = " << a << endl;
	//cout << "n = " << n << endl;
	//Calculamos a = a (mod p)
	mpz_mod(a2.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
	//cout << "a2 = " << a2 << endl << flush;
	
	//Comprobaciones iniciales
    if (a2 == 0) {
        sqrt_a = 0;
    }
    else{
    	if (n == 2) {
        	sqrt_a = a2 % n;
    	}
    	else{
			if (mpz_legendre(a2.get_mpz_t(), n.get_mpz_t()) != 1) {
		 	   throw std::invalid_argument("No existe una raíz cuadrada modular");
			}
			else{
				//Calculamos q,s tales que a = q*2^s, con q impar
				while(q % 2 == 0){
					q = q / 2;
					s++;
				}
				//Buscamos un numero z que no sea residuo cuadratico en Z_n
				z = 2;
				while(mpz_legendre(z.get_mpz_t(),n.get_mpz_t()) != -1){
					z++;
				}
				//cout << "S = " << s << endl;
				//cout << "Q = " << q << endl;
				
				//cout << "Valores iniciales" << endl;
				M = s;
				mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
				mpz_powm(t.get_mpz_t(), a2.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
				aux = (q+1)/2;
				mpz_powm(R.get_mpz_t(), a2.get_mpz_t(), aux.get_mpz_t(), n.get_mpz_t());
				/*
				cout << "M = " << M << endl;
				cout << "c = " << c << endl;
				cout << "aux = " << aux << endl;
				cout << "t = " << t << endl;
				cout << "R = " << R << endl;
				*/
				while (t != 0 && t != 1) {
					//cout << "Iteracion" << endl;
					mpz_class tt = t;
					unsigned int i = 0;
					for (i = 0; i < M.get_ui(); i++) {
						if (tt == 1) {
						    break;
						}
						mpz_powm_ui(tt.get_mpz_t(), tt.get_mpz_t(), 2, n.get_mpz_t());
					}

					if (i == M.get_ui()) {
						throw std::invalid_argument("No existe una raíz cuadrada modular");
					}

					mpz_powm_ui(b.get_mpz_t(), c.get_mpz_t(), 1 << (M.get_ui() - i - 1), n.get_mpz_t());
					//cout << "b = " << b << endl;
					
					M = i;
					//cout << "M = " << M << endl;
					
					//c = b^2 (mod n)
					c = (b*b) % n;
					//cout << "c = " << c << endl;
					
					//R = R*b (mod n)
					R = (R*b) % n;
					//cout << "R = " << R << endl;
					
					//t = t*b^2 (mod n) = t*c (mod n)
					t = (t*c) % n;
					//cout << "t = " << t << endl;
				}
				if(t == 0){
					sqrt_a = 0;
				}
				else{
					sqrt_a = R;
				}
				
			}
		}
	}
}

// Función para transponer una matriz
vector<vector<bool>> transpose(vector<vector<bool>>& matrix) {
    int n = matrix.size();
    int m = matrix[0].size();
    vector<vector<bool>> transposed(m, vector<bool>(n));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}


// Función para realizar la eliminación gaussiana sobre matrix2 en mod 2
void gaussian_elimination(vector<vector<bool>>& matrix) {
    int n = matrix.size();
    int m = matrix[0].size();
    
    int row = 0;
    for (int col = 0; col < m && row < n; ++col) {
        // Buscamos una fila con un 1 en la columna col
        int sel = row;
        for (int i = row; i < n; ++i) {
            if (matrix[i][col]) {
                sel = i;
                break;
            }
        }
        
        // Si no hay ningún 1 en la columna col, seguimos
        if (!matrix[sel][col]) continue;

        // Intercambiamos la fila seleccionada con la fila actual
        swap(matrix[sel], matrix[row]);

        // Escribimos ceros en todas las filas excepto la actual
        for (int i = 0; i < n; ++i) {
            if (i != row && matrix[i][col]) {
                for (int j = col; j < m; ++j) {
                    matrix[i][j] = matrix[i][j] ^ matrix[row][j];
                }
            }
        }
        row++;
    }
}

// Funcion para encontrar las soluciones del sistema anterior
vector<vector<bool>> find_solutions(const vector<vector<bool>>& matrix) {
    int n = matrix.size();
    int m = matrix[0].size();
    
    vector<int> pivot(m, -1); // Índice de fila del pivote para cada columna, -1 si no hay pivote
    vector<vector<bool>> solutions;

    // Identificamos columnas pivote
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (matrix[i][j]) {
                pivot[j] = i;
                break;
            }
        }
    }

    // Generamos todas las combinaciones de variables libres
    int num_free_vars = count(pivot.begin(), pivot.end(), -1);
    int num_solutions = 1 << num_free_vars; // 2^num_free_vars combinaciones posibles

    for (int k = 0; k < num_solutions; ++k) {
        vector<bool> solution(m, false);
        int free_var_idx = 0;

        // Asignar valores a variables libres según la combinación k
        for (int j = 0; j < m; ++j) {
            if (pivot[j] == -1) {
                solution[j] = (k >> free_var_idx) & 1;
                free_var_idx++;
            }
        }

        // Calculamos los valores de las variables pivote
        for (int j = 0; j < m; ++j) {
            if (pivot[j] != -1) {
                bool value = false;
                for (int l = j + 1; l < m; ++l) {
                    if (matrix[pivot[j]][l]) {
                        value ^= solution[l];
                    }
                }
                solution[j] = value;
            }
        }

        solutions.push_back(solution);
    }

    return solutions;
}


mpz_class factorizacionCribaCuadratica(const mpz_class &n, gmp_randstate_t state, const mpz_class &k, const unsigned int &tam_tabla){
	vector<mpz_class> tabla_criba;
	vector<mpz_class> base_primos;
	vector<mpz_class> x_uniformes;
	vector<vector<bool>> factores_uniformes;
	vector<vector<int>> factores_uniformes_count;
	mpz_class p_actual = 3;
	mpz_class x_ini, entrada, r1, r2, sol1, sol2, it, p, q;
	int s_legendre;
	bool exito = false;
	
	//Buscamos nuestra base de primos
	//Incluiremos en ella aquellos primos p: n % p = y^2 (es decir, n es resto cuadratico modulo p)
	//Esto es equivalente a ver si el simbolo de Legendre de n/p es igual a 1.
	base_primos.push_back(2);
	while(k > p_actual){
		s_legendre = mpz_legendre(n.get_mpz_t(),p_actual.get_mpz_t());
		if(s_legendre == 1){
			base_primos.push_back(p_actual);
			//cout << p_actual << endl;
		}
		mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
	}
	for(int i = 0; i < base_primos.size(); i++){
		//cout << base_primos[i] << endl;
	}
	
	
	//Creamos la tabla de la criba y la tabla de factores
	mpz_sqrt(x_ini.get_mpz_t(), n.get_mpz_t());	//Guardamos en x_ini = sqrt(n)
	x_ini++;
	for(int i = 0; i < tam_tabla; i++){
		entrada = i + x_ini;
		entrada = entrada*entrada - n;
		tabla_criba.push_back(entrada);
		//cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << endl;
	}
	
	vector<vector<bool>> factores(tam_tabla, vector<bool>(base_primos.size(), false));
	vector<vector<int>> factores_count(tam_tabla, vector<int>(base_primos.size(), 0));
	
	//Calculamos las soluciones de f(x) = x^2 - n = 0 (mod p)
	//--> x^2 = n (mod p) --> x = sqrt(n) (mod p)
	//Una vez las tengamos, ejecutamos la criba
	for(int i = 0; i < base_primos.size(); i++){
		//Calculamos la raiz cuadrada de n
		//cout << "Primo " << base_primos[i] << " " << endl;
		sqrt_mod(r1, n, base_primos[i]);
		
		
		//Realizamos la criba con la solucion positiva
		sol1 = (r1 - x_ini);
		mpz_mod(sol1.get_mpz_t(), sol1.get_mpz_t(), base_primos[i].get_mpz_t());
		//cout << "criba con " << sol1;
		it = sol1;
		while(it < tabla_criba.size()){
			tabla_criba[it.get_ui()] /= base_primos[i];
			factores[it.get_ui()][i] = !factores[it.get_ui()][i];
			factores_count[it.get_ui()][i] += 1;
			it += base_primos[i];
		}
		
		//Calculamos la solucion negativa
		r2 = (-r1) % base_primos[i];
		sol2 = (r2 - x_ini);
		mpz_mod(sol2.get_mpz_t(), sol2.get_mpz_t(), base_primos[i].get_mpz_t());
		//Realizamos la criba con la solucion negativa
		if(sol1 != sol2){
			//cout << " y criba con " << sol2;
			while(sol2 < tabla_criba.size()){
				tabla_criba[sol2.get_ui()] /= base_primos[i];
				factores[sol2.get_ui()][i] = !factores[sol2.get_ui()][i];
				factores_count[sol2.get_ui()][i] += 1;
				sol2 += base_primos[i];
			}
		}
		//cout << endl;
	}
	
	//Extraemos los restos uniformes (=1)
	for(int i = 0; i < tabla_criba.size(); i++){
		if(tabla_criba[i] == 1){
			x_uniformes.push_back(i);
			factores_uniformes.push_back(factores[i]);
			factores_uniformes_count.push_back(factores_count[i]);
		}
		/*
		cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << "\t[";
		
		for(int j = 0; j < base_primos.size(); j++){
			//cout << factores[i][j] << "\t";
			cout << factores_count[i][j] << "\t";
		}
		cout << "]" << endl;
		*/
	}	
	
	for(int i = 0; i < x_uniformes.size(); i++){
		cout << "Resto " << i << " = " << x_uniformes[i] << "\t[";
		for(int j = 0; j < base_primos.size() - 1; j++){
			//cout << factores_uniformes[i][j] << "\t";
			cout << factores_uniformes_count[i][j] << "\t";
		}
		cout << factores_uniformes[i][base_primos.size() - 1] << "]" << "\t";
		cout << factores_uniformes_count[i][base_primos.size() - 1] << endl;
	}
	
	//Trasponemos factores_uniformes y ejecutamos la eliminacion gaussiana
	vector<vector<bool>> factores_uniformes_t = transpose(factores_uniformes);
	
	gaussian_elimination(factores_uniformes_t);
	//cout << "Eliminacion gaussiana" << endl;
	/*
	for(int i = 0; i < factores_uniformes_t.size(); i++){
		for(int j = 0; j < factores_uniformes_t[0].size(); j++){
			//cout << factores_uniformes_t[i][j] << "\t";
		}
		//cout << endl;
	}
	*/
	vector<vector<bool>> sols = find_solutions(factores_uniformes_t);
	cout << sols.size() << endl;
	vector<bool> sol;
	if(sols.size() > 0){
		for(int n_sol = 0; n_sol < sols.size() && !exito; n_sol++){
			sol = sols[n_sol];
			cout << "Solucion \n[";
			for(int i = 0; i < sol.size(); i++){
				cout << sol[i] << "\t";
			}
			cout << "]" << endl;
			
			vector<int> factores_r(base_primos.size(), false);
			for(int i = 0; i < sol.size(); i++){
				if(sol[i]){
					for(int j = 0; j < base_primos.size(); j++){
						factores_r[j] += factores_uniformes_count[i][j];
					}
				}
			}
			
			
			cout << "COMPROBACION" << endl << "\t";
			for(int i = 0; i < factores_r.size(); i++){
				cout << factores_r[i] << "\t";
			}
			cout << endl;
			
			for(int j = 0; j < base_primos.size(); j++){
				factores_r[j] /= 2;
			}
			
			//Calculamos la congruencia resultante x^2 = y^2 mod n
			mpz_class x_acum = 1;
			mpz_class r_acum = 1, r2_acum = 1;
			mpz_class x_actual;
			for(int i = 0; i < sol.size(); i++){
				//cout << "i = " << i << endl;
				if(sol[i] == 1){
					x_actual = x_ini + x_uniformes[i];
					//cout << "x_uniformes[i] = " << x_uniformes[i] << endl;
					cout << x_uniformes[i] << "\t";
					//cout << "Resto " << i << " = " << x_uniformes[i] << "\t[";
					for(int j = 0; j < base_primos.size() - 1; j++){
						//cout << factores_uniformes[i][j] << "\t";
						cout << factores_uniformes_count[i][j] << "\t";
					}
					//cout << factores_uniformes[i][base_primos.size() - 1] << "]" << "\t";
					cout << factores_uniformes_count[i][base_primos.size() - 1] << endl;
					
					x_acum *= x_actual;
				}
			}
			
			for(int j = 0; j < factores_r.size(); j++){
				for(int k = 0; k < factores_r[j]; k++)
				r_acum *= base_primos[j];
			}
			
			cout << "Congruencia " << endl;
			
			mpz_mod(x_acum.get_mpz_t(), x_acum.get_mpz_t(), n.get_mpz_t());
			mpz_mod(r_acum.get_mpz_t(), r_acum.get_mpz_t(), n.get_mpz_t());
			cout << "x = " << x_acum << endl;
			cout << "x^2 = " << x_acum*x_acum % n << endl;
			cout << "y = " << r_acum << endl;
			cout << "y^2 = " << r_acum*r_acum % n << endl;
			
			//Calculamos el mcd de la diferencia de ambos
			mpz_class dif = (x_acum - r_acum) % n;
			mpz_gcd(p.get_mpz_t(), dif.get_mpz_t(), n.get_mpz_t());
			q = n / p;
			cout << "p = " << p << endl;
			cout << "q = " << q << endl;
			if(p != 1 && q != 1)
				exito = true;
		}
	}
	else{
		cout << "El sistema no tiene solucion :(" << endl;
	}
	return p;
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
	
	bool exito;
	unsigned int att = 1000;
	mpz_class n, p, q, k;
	
	//mpz_set_ui(n, 265905204186593);
	n = 265905204186593;
	k = 30;
	p = factorizacionCurvasElipticas(n, state, k, att);
	if(p != n){
		cout << "n = " << n << endl;
		cout << "p = " << p << endl;
		cout << "q = " << q << endl;
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
	p = factorizacionCribaCuadratica(n1, state, k1, 1000000);
	
    return 0;
}

