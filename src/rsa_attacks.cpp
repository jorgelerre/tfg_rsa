#include "utils.h"
#include "rsa_attacks.h"
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

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
mpz_class rhoPollard(mpz_class n, gmp_randstate_t state, bool debug){
	mpz_class p, t, l, dif, mcd = 1;
	
	//Comprobamos que no sea primo
	if(millerRabin(n, 10, state)){
		cout << "p es primo" << endl;
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
			if(debug){
				cout << "t = " << t << endl;
				cout << "l = " << l << endl;
				cout << "|t - l| = " << dif << endl << flush;
			}
			mpz_gcd (mcd.get_mpz_t(), n.get_mpz_t(), dif.get_mpz_t());	//Calculamos el mcd
		}
		p = mcd;
		if(debug) cout << "p = " << p << endl; 
	}
	
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

mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k,
									   const unsigned int att, bool debug){
	Punto Q, M;
	mpz_class p, a, b, aux, inv, L;
	bool exito = false;
	
	//Comprobamos que n no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		for(unsigned int i = 0; i < att && !exito; i++){
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
			
			if(debug){
				cout << "Intento " << att << endl;
				cout << "Curva eliptica: " << endl;
				cout << "a = " << a << endl;
				cout << "b = " << b << endl;
				cout << "Punto en la curva = (" << Q.x << "," << Q.y << ")" << endl;
				cout << k << "-potencia-suave = " << L << endl;
			}
			
			//Realizamos la multiplicacion
			M = multiplicacionCurvaEliptica(Q, L, a, b, n, inv);
			
			//Si hemos encontrado un elemento no invertible, obtenemos un factor de n
			if(M.x == -2){
				
				//Calculamos gcd
				mpz_gcd (p.get_mpz_t(), inv.get_mpz_t(), n.get_mpz_t());
				exito = true;
				
				if(debug){
					cout << "Encontramos un elemento no invertible!!" << endl;
					cout << "inv = " << inv << endl;
					cout << "p = " << p << endl;
				}
			}
			else{
				if(debug){
					if(M.x == -1){
						cout << "Encontramos una solucion trivial." << endl;
					}
					else{
						cout << "No hemos encontrado solucion." << endl;
					}
				}
			}
		}
	}
	//Si no hemos conseguido factorizar, devolvemos p = n
	if(!exito){
		if(debug) cout << "Iteraciones agotadas, devolviendo solucion trivial (p = n)" << endl;
		p = n;
	}
	
	return p;
}


// Función para encontrar una raíz cuadrada de a módulo p usando el algoritmo de Tonelli-Shanks
mpz_class sqrt_mod(const mpz_class &a, const mpz_class &n, bool debug){
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
vector<vector<bool>> gaussian_elimination(vector<vector<bool>> &matrix) {
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


mpz_class factorizacionCribaCuadratica(const mpz_class &n, const mpz_class &k, 
									   const unsigned int tam_tabla, bool debug){
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
			if(debug) cout << p_actual << endl;
		}
		mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
	}
	for(unsigned int i = 0; i < base_primos.size(); i++){
		if(debug) cout << base_primos[i] << endl;
	}
	
	
	//Creamos la tabla de la criba y la tabla de factores
	mpz_sqrt(x_ini.get_mpz_t(), n.get_mpz_t());	//Guardamos en x_ini = sqrt(n)
	x_ini++;
	for(unsigned int i = 0; i < tam_tabla; i++){
		entrada = i + x_ini;
		entrada = entrada*entrada - n;
		tabla_criba.push_back(entrada);
		if(debug) cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << endl;
	}
	vector<vector<bool>> factores(tam_tabla, vector<bool>(base_primos.size(), false));
	vector<vector<int>> factores_count(tam_tabla, vector<int>(base_primos.size(), 0));
	//Calculamos las soluciones de f(x) = x^2 - n = 0 (mod p)
	//--> x^2 = n (mod p) --> x = sqrt(n) (mod p)
	//Una vez las tengamos, ejecutamos la criba
	for(unsigned int i = 0; i < base_primos.size(); i++){
		//Calculamos la raiz cuadrada de n
		if(debug) cout << "Primo " << base_primos[i] << " " << endl;
		r1 = sqrt_mod(n, base_primos[i]);
		
		
		//Realizamos la criba con la solucion positiva
		sol1 = (r1 - x_ini);
		mpz_mod(sol1.get_mpz_t(), sol1.get_mpz_t(), base_primos[i].get_mpz_t());
		if(debug) cout << "criba con " << sol1;
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
			if(debug) cout << " y criba con " << sol2;
			while(sol2 < tabla_criba.size()){
				tabla_criba[sol2.get_ui()] /= base_primos[i];
				factores[sol2.get_ui()][i] = !factores[sol2.get_ui()][i];
				factores_count[sol2.get_ui()][i] += 1;
				sol2 += base_primos[i];
			}
		}
		if(debug) cout << endl;
	}
	
	
	
	for(unsigned int i = 0; i < tabla_criba.size(); i++){
		if(tabla_criba[i] == 1){
			x_uniformes.push_back(i);
			factores_uniformes.push_back(factores[i]);
			factores_uniformes_count.push_back(factores_count[i]);
		}
			if(debug) {
			cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << "\t[";
			
			for(unsigned int j = 0; j < base_primos.size(); j++){
				//cout << factores[i][j] << "\t";
				cout << factores_count[i][j] << "\t";
			}
			cout << "]" << endl;
			//Extraemos los restos uniformes (=1)
		}
	}	
	cout << "NEXT" << x_uniformes.size() << endl;
	if(x_uniformes.size() > 0){
		if(debug) {
			for(unsigned int i = 0; i < x_uniformes.size(); i++){
				cout << "Resto " << i << " = " << x_uniformes[i] << "\t[";
				for(unsigned int j = 0; j < base_primos.size() - 1; j++){
					//cout << factores_uniformes[i][j] << "\t";
					cout << factores_uniformes_count[i][j] << "\t";
				}
				cout << factores_uniformes[i][base_primos.size() - 1] << "]" << "\t";
				cout << factores_uniformes_count[i][base_primos.size() - 1] << endl;
			}
		}
		//Trasponemos factores_uniformes y ejecutamos la eliminacion gaussiana
		vector<vector<bool>> factores_uniformes_t = transpose(factores_uniformes);
		
		factores_uniformes_t = gaussian_elimination(factores_uniformes_t);
		if(debug){
			cout << "Eliminacion gaussiana" << endl;
			for(unsigned int i = 0; i < factores_uniformes_t.size(); i++){
				for(unsigned int j = 0; j < factores_uniformes_t[0].size(); j++){
					cout << factores_uniformes_t[i][j] << "\t";
				}
				cout << endl;
			}
		}
		vector<vector<bool>> sols = find_solutions(factores_uniformes_t);
		if(debug) cout << sols.size() << endl;
		vector<bool> sol;
		if(sols.size() > 0){
			for(unsigned int n_sol = 0; n_sol < sols.size() && !exito; n_sol++){
				sol = sols[n_sol];
				if(debug) {
					cout << "Solucion \n[";
					for(unsigned int i = 0; i < sol.size(); i++){
						cout << sol[i] << "\t";
					}
					cout << "]" << endl;
				}
				vector<int> factores_r(base_primos.size(), false);
				for(unsigned int i = 0; i < sol.size(); i++){
					if(sol[i]){
						for(unsigned int j = 0; j < base_primos.size(); j++){
							factores_r[j] += factores_uniformes_count[i][j];
						}
					}
				}
				
				if(debug){
					cout << "COMPROBACION" << endl << "\t";
					for(unsigned int i = 0; i < factores_r.size(); i++){
						cout << factores_r[i] << "\t";
					}
					cout << endl;
				}
				for(unsigned int j = 0; j < base_primos.size(); j++){
					factores_r[j] /= 2;
				}
				
				//Calculamos la congruencia resultante x^2 = y^2 mod n
				mpz_class x_acum = 1;
				mpz_class r_acum = 1, r2_acum = 1;
				mpz_class x_actual;
				for(unsigned int i = 0; i < sol.size(); i++){
					if(debug) cout << "i = " << i << endl;
					if(sol[i] == 1){
						x_actual = x_ini + x_uniformes[i];
						if(debug) {
							cout << "x_uniformes[i] = " << x_uniformes[i] << endl;
							cout << x_uniformes[i] << "\t";
							//cout << "Resto " << i << " = " << x_uniformes[i] << "\t[";
							
							for(unsigned int j = 0; j < base_primos.size() - 1; j++){
								//cout << factores_uniformes[i][j] << "\t";
								cout << factores_uniformes_count[i][j] << "\t";
							}
							
							//cout << factores_uniformes[i][base_primos.size() - 1] << "]" << "\t";
							cout << factores_uniformes_count[i][base_primos.size() - 1] << endl;
						}
						x_acum *= x_actual;
					}
				}
				
				for(unsigned int j = 0; j < factores_r.size(); j++){
					for(int k = 0; k < factores_r[j]; k++)
					r_acum *= base_primos[j];
				}
				
				if(debug) cout << "Congruencia " << endl;
				
				mpz_mod(x_acum.get_mpz_t(), x_acum.get_mpz_t(), n.get_mpz_t());
				mpz_mod(r_acum.get_mpz_t(), r_acum.get_mpz_t(), n.get_mpz_t());
				if(debug) {
					cout << "x = " << x_acum << endl;
					cout << "x^2 = " << x_acum*x_acum % n << endl;
					cout << "y = " << r_acum << endl;
					cout << "y^2 = " << r_acum*r_acum % n << endl;
				}
				//Calculamos el mcd de la diferencia de ambos
				mpz_class dif = (x_acum - r_acum) % n;
				mpz_gcd(p.get_mpz_t(), dif.get_mpz_t(), n.get_mpz_t());
				q = n / p;
				if(debug){
					cout << "p = " << p << endl;
					cout << "q = " << q << endl;
				}
				if(p != 1 && q != 1)
					exito = true;
			}
		}
		else{
			if(debug) cout << "El sistema no tiene solucion :(" << endl;
			p = n;
		}
	}
	else{
		if(debug) cout << "No se ha encontrado ningun resto k-uniforme :(" << endl;
		p = n;
	}
	return p;
}
