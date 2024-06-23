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
unsigned long int generateRandomLimb(gmp_randstate_t state) {
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
	bool no_fail = true, primo_bajo = false;
	mpz_class b;
	
	//Realizamos el test de las divisiones sucesivas con los primeros primos (<8 bits)
	const int NUM_PRIMES = 54;
	const unsigned int primes[NUM_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 
			   59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 
	 		   157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};
	
	for(int i = 0; i < NUM_PRIMES && no_fail; i++){
		b = p % primes[i];
		if(b == 0){
			if(p == primes[i]){	//Si p es igual al primo, devolvemos que es primo
				no_fail = true;
				primo_bajo = true;
			}	
			else{	//Si no, p es multiplo del primo, y por tanto es compuesto
				no_fail = false;
			}
		}	
	}
	
	//Realizamos it veces la comprobacion de Miller-Rabin
	if(!primo_bajo){
		for(int i = 0; i < it && no_fail; i++){
			//Generamos un numero aleatorio (pequeño (1 miembro) para aligerar los calculos)
			b = generateRandomLimb(state);
			//Comprobamos la pseudoprimalidad en dicha base
			no_fail = strongPseudoprime(p, b);
		}
	}
	return no_fail;
}

/**
 * Función para generar un número primo de 'bits' bits
 * @bits Numero de bits deseado del primo a generar
 * @state Estado del generador de numeros aleatorios de GMP
 * @return Un primo del tamanio especificado
 */
mpz_class generatePrime(unsigned int bits, gmp_randstate_t state) {
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
mpz_class generateStrongPrime(unsigned int bits, gmp_randstate_t state, bool debug) {
    mpz_class p,p_0,r,s,t,i,j,aux;
    int bits_j, bits_1, bits_2;
    
    if(debug) cout << "Buscando s,t..." << endl;
    //Generamos dos primos grandes s, t
    if(bits > 90){
		bits_1 = (bits-log2(bits))/2 - 3;
		bits_2 = bits_1 - log2(bits_1) - 6;
    }
    else{
    	bits_1 = (bits-log2(bits))/2 - 1;
		bits_2 = bits_1 - log2(bits_1) - 1;
		cout << "bits 1: " << bits_1 << endl;
		cout << "bits 2: " << bits_2 << endl;
		if(bits_1 < 4){
			bits_1 = 4;
		}
		if(bits_2 < 2){
			bits_2 = 2;
		}
		
		cout << "bits 1: " << bits_1 << endl;
		cout << "bits 2: " << bits_2 << endl;
    }
    
    do{
		s = generatePrime(bits_1, state);
		t = generatePrime(bits_2, state);
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
    if(debug) cout << "primo final = " << p << endl;
    return p;
}



//Suma en la curva eliptica y^2 = x^3 + ax + b en el cuerpo Z_p
Punto sumaCurvaEliptica(const Punto s1, const Punto s2, const mpz_class a, 
					    const mpz_class b, const mpz_class p, mpz_class &inv, bool debug){
	Punto s3;
	mpz_class aux, lambda;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente
	if(debug){
		cout << "Suma curva eliptica" << endl;
		cout << "s1 = (" <<  s1.x << ", " 
						  <<  s1.y << ")" << endl;		  
		cout << "s2 = (" << s2.x << ", " 
						  << s2.y << ")" << endl;
		cout << "p = " << p << endl;
	}
	//Si alguno de los puntos es O, devolvemos el otro punto como solucion
	if(s1.x == 0){
		s3.x = s2.x;
		s3.y = s2.y;
		if(debug) cout << "Sumando 1 es \"O\": devolviendo s2" << endl;
	}
	else{
		if(s2.x < 0){
			s3.x = s1.x;
			s3.y = s1.y;
			if(debug) cout << "Sumando 2 es \"O\": devolviendo s1" << endl;
		}
		else{
			//Si tenemos dos puntos de la forma (x,y) y (x,-y), devolvemos O
			aux = s1.y + s2.y;
			if(debug) cout << "Suma curva eliptica sin O" << endl;
			mpz_mod(aux.get_mpz_t(), aux.get_mpz_t(), p.get_mpz_t());
			if(debug) cout << "Suma curva eliptica sin O" << endl;
			if(s1.x == s2.x && aux == 0){
				s3.x = -1;
				s3.y = -1;
				if(debug) cout << "Puntos opuestos: devolviendo \"O\"" << endl;
			}
			//En otro caso
			else{
				//Calculamos lambda
				//Si ambos puntos son iguales
				if(s1.x == s2.x && s1.y == s2.y){
					if(debug) cout << "Puntos iguales: calculando lambda" << endl;
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
					if(debug) cout << "Puntos diferentes: calculando lambda" << endl;
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
					if(debug) cout << "lambda = " << lambda << endl;
					//Calculamos la coordenada x de la suma
					//x3 = lambda^2 - x1 - x2 (mod p)
					s3.x = lambda*lambda - s1.x - s2.x;
					mpz_mod(s3.x.get_mpz_t(), s3.x.get_mpz_t(), p.get_mpz_t());
					
					//Calculamos la coordenada y
					//y3 = lambda*(x_1-x_3) - y1 (mod p)
					s3.y = lambda*(s1.x - s3.x) - s1.y;
					mpz_mod(s3.y.get_mpz_t(), s3.y.get_mpz_t(), p.get_mpz_t());
					if(debug) cout << "Suma = (" << s3.x << "," << s3.y << ")" << endl;
				}
				else{
					//Punto no definido
					if(debug){
						cout << "Suma no definida: no existe el inverso de " << inv << endl;
						cout << "inv = " << inv << endl;
					}
					s3.x = -2;
					s3.y = -2;
				}
			}
		}
	}
	return s3;
}

Punto multiplicacionCurvaEliptica(const Punto f1, const mpz_class f2, const mpz_class a, 
								const mpz_class b, const mpz_class p, mpz_class &inv, bool debug){
	Punto p2, mul;
	mpz_class r, lambda, f2_2;
	int exito = 1;	//Comprueba si el calculo de la inversa se hizo correctamente en las sumas
	int i = 0;
	
	if(debug) cout << "Inicio multiplicacion" << endl;
	
	//Inicializamos la multiplicacion en "O"
	mul.x = -1;
	mul.y = -1;
	
	//Guardamos el punto a multiplicar P en p2
	p2.x = f1.x;
	p2.y = f1.y;
	f2_2 = f2;
	if(debug) cout << "p2 = (" << p2.x << "," << p2.y << ")" << endl;
	
	while(f2_2 > 0 && exito != 0){
		if(debug){
			cout << "Iteracion " << i << endl;
			cout << "f2 = " << f2_2 << endl;
		}
		//Calculamos el siguiente bit
		mpz_fdiv_qr_ui(f2_2.get_mpz_t(), r.get_mpz_t(), f2_2.get_mpz_t(), 2);
		if(debug) cout << "Valor bit actual: " << r << endl;
		//Si el bit es 1, multiplicamos p2 por mul
		if(r == 1){
			mul = sumaCurvaEliptica(mul, p2, a, b, p, inv);
			if(mul.x == -2){
				exito = 0;
			}
		}
		//Calculamos la siguiente potencia de 2 de P
		if(exito != 0){
			p2 = sumaCurvaEliptica(p2, p2, a, b, p, inv);
			if(p2.x == -2){
				exito = 0;
			}
		}
		if(debug) cout << "mul actual = (" << mul.x << "," << mul.y << ")" << endl;
		i++;
	}
	//Devolvemos un punto no existente si falla alguna suma
	if(exito == 0){
		if(debug) cout << "La multiplicacion no queda definida: no existe el inverso de " << inv << endl;
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

// Función para encontrar una raíz cuadrada de a módulo p usando el algoritmo de Tonelli-Shanks
mpz_class sqrtMod(const mpz_class &a, const mpz_class &n, bool debug){
	mpz_class sqrt_a, a2, b, q = n-1, s = 0, z, aux;
	mpz_class M, c, t, R;
	
	//Calculamos a = a (mod p)
	mpz_mod(a2.get_mpz_t(), a.get_mpz_t(), n.get_mpz_t());
	if(debug){
		cout << "a = " << a << endl;
		cout << "n = " << n << endl;
		cout << "a2 = " << a2 << endl << flush;
	}
	//Comprobaciones iniciales
    if (a2 == 0) {
    	if(debug) cout << "La raiz de 0 siempre es 0" << endl;
        sqrt_a = 0;
    }
    else{
    	if (n == 2) {
    		if(debug) cout << "La raiz de en Z_2 siempre es 0 sii es par" << endl;
        	sqrt_a = a2 % n;
    	}
    	else{
			if (mpz_legendre(a2.get_mpz_t(), n.get_mpz_t()) != 1) {
				if(debug){
					cout << "No existe la raiz cuadrada modular: el simbolo de Legendre no es 1." << endl;
				}
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
				if(debug){
					cout << "S = " << s << endl;
					cout << "Q = " << q << endl;
					cout << "Valores iniciales" << endl;
				}
				M = s;
				mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
				mpz_powm(t.get_mpz_t(), a2.get_mpz_t(), q.get_mpz_t(), n.get_mpz_t());
				aux = (q+1)/2;
				mpz_powm(R.get_mpz_t(), a2.get_mpz_t(), aux.get_mpz_t(), n.get_mpz_t());
				if(debug){
					cout << "M = " << M << endl;
					cout << "c = " << c << endl;
					cout << "aux = " << aux << endl;
					cout << "t = " << t << endl;
					cout << "R = " << R << endl;
				}
				while (t != 0 && t != 1) {
					if(debug) cout << "Iteracion" << endl;
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
					M = i;
					c = (b*b) % n;
					R = (R*b) % n;
					//t = t*b^2 (mod n) = t*c (mod n)
					t = (t*c) % n;
					if(debug){
						cout << "-----------------" << endl;
						cout << "b = " << b << endl;
						cout << "M = " << M << endl;
						cout << "c = " << c << endl;
						cout << "R = " << R << endl;
						cout << "t = " << t << endl;
					}
					
				}
				if(t == 0){
					sqrt_a = 0;
				}
				else{
					sqrt_a = R;
				}
				if(debug) cout << "Raiz cuadrada = " << sqrt_a << endl;
			}
		}
	}
	return sqrt_a;
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
vector<vector<bool>> eliminacionGaussiana(const vector<vector<bool>> &matrix) {
    vector<vector<bool>> matrix2 = matrix;
    int n = matrix2.size();
    int m = matrix2[0].size();
    
    int row = 0;
    for (int col = 0; col < m && row < n; ++col) {
        // Buscamos una fila con un 1 en la columna col
        int sel = row;
        for (int i = row; i < n; ++i) {
            if (matrix2[i][col]) {
                sel = i;
                break;
            }
        }
        
        // Si no hay ningún 1 en la columna col, seguimos
        if (!matrix2[sel][col]) continue;

        // Intercambiamos la fila seleccionada con la fila actual
        swap(matrix2[sel], matrix2[row]);

        // Escribimos ceros en todas las filas excepto la actual
        for (int i = 0; i < n; ++i) {
            if (i != row && matrix2[i][col]) {
                for (int j = col; j < m; ++j) {
                    matrix2[i][j] = matrix2[i][j] ^ matrix2[row][j];
                }
            }
        }
        row++;
    }
    return matrix2;
}

// Funcion para encontrar las soluciones del sistema anterior
vector<vector<bool>> encuentraSoluciones(const vector<vector<bool>>& matrix) {
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
	
	if(num_solutions > 100){	//Si hay mas de 100 soluciones, nos quedamos con las primeras solo
		num_solutions = 100;
	}
	
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

//Calcula la fraccion continua de la fraccion num/den
vector<mpz_class> cocientesFraccionContinua(const mpz_class &num, const mpz_class &den) {
    mpz_class n = num, d = den, q, r;
    vector<mpz_class> cf;
    while (d != 0) {
    	mpz_fdiv_qr(q.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t(), d.get_mpz_t());
        cf.push_back(q);
        n = d;
        d = r;
    }
    return cf;
}

//Resuelve ecuacion cuadratica con a y c no nulos
bool resuelveEcuacionCuadratica(mpz_class &p, mpz_class &q, mpz_class a, mpz_class b, mpz_class c){
	mpz_class discriminant = b * b - 4 * a * c;
	bool exito = false;
	if (discriminant >= 0) {
		mpz_class sqrt_discriminant;
		mpz_sqrt(sqrt_discriminant.get_mpz_t(), discriminant.get_mpz_t());
		p = (-b + sqrt_discriminant) / (2 * a);
		q = (-b - sqrt_discriminant) / (2 * a);
		exito = true;
	}
	return exito;
}
