#include "utils.h"
#include "rsa_attacks.h"
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;


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

mpz_class p1Pollard(mpz_class n, gmp_randstate_t state, const mpz_class k, const unsigned int att, bool debug){
	mpz_class a, B, x, mcd = 1, p_potencia, p_actual;
	bool exito = false;
	
	for(unsigned int i = 0; i < att && !exito; i++){
		//Seleccionamos una base a
		mpz_urandomm (a.get_mpz_t(), state, n.get_mpz_t());
		
		p_actual = 2;
		//Mientras el primo actual no supere el valor de K
		while(p_actual < k && !exito){
			cout << p_actual << endl;
			p_potencia = p_actual;
			//Calculamos el valor de p_actual^e tal que no supere K
			while(p_potencia < k){
				p_potencia = p_potencia * p_actual;	//p_potencia = p_actual^{n+1}
			}
			//Dividimos para que p_potencia < k
			mpz_divexact(p_potencia.get_mpz_t(), p_potencia.get_mpz_t(), p_actual.get_mpz_t());	
			//Multiplicamos dicho valor en el acumulador
			mpz_powm(a.get_mpz_t(), a.get_mpz_t(), p_potencia.get_mpz_t(), n.get_mpz_t());
			//Calculamos el mcd de a^B-1 y n
			x = a - 1;
			mpz_gcd(mcd.get_mpz_t(), x.get_mpz_t(), n.get_mpz_t());
			
			//Si mcd no es trivial, hemos conseguido un factor de n
			if(mcd < n && mcd > 1)
				exito = true;
			//Pasamos al siguiente primo
			mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
		}
		
		if(debug){
			if(mcd == n){
				cout << "mcd = n : k es demasiado alto." << endl;
			}
			else{
				if(mcd == 1){
					cout << "mcd = 1 : k es demasiado bajo o la base no es la adecuada." << endl;
				}
				else{
					cout << "mcd = " << mcd << ". Hemos factorizado n." << endl;
				}
			}
		}
	}
	
	return mcd;
}


mpz_class factorizacionCurvasElipticas(const mpz_class n, gmp_randstate_t state, const mpz_class k,
									   const unsigned int att, bool debug){
	Punto Q, M;
	mpz_class p, a, b, aux, inv, p_actual, p_potencia;
	bool exito = false;
			
	//Comprobamos que n no sea primo
	if(millerRabin(n, 10, state)){
		p = n;
	}
	else{
		for(unsigned int i = 0; i < att && !exito; i++){
			//Escogemos una pseudo curva eliptica eligiendo a, b y definiendola en Z_n
			mpz_urandomm (a.get_mpz_t(), state, n.get_mpz_t());
			
			//Escogemos un punto Q aleatorio de la curva
			mpz_urandomm (Q.x.get_mpz_t(), state, n.get_mpz_t());	//Elegimos x aleatorio
			mpz_urandomm (Q.y.get_mpz_t(), state, n.get_mpz_t());
			//Calculamos b como y^2-x^3 - ax (mod n)
			b = Q.y * Q.y - Q.x * Q.x * Q.x + a*Q.x;
			
			mpz_mod(b.get_mpz_t(), b.get_mpz_t(), n.get_mpz_t());
			
			if(debug){
				cout << "Intento " << i << endl;
				cout << "Curva eliptica: " << endl;
				cout << "a = " << a << endl;
				cout << "b = " << b << endl;
				cout << "Punto en la curva = (" << Q.x << "," << Q.y << ")" << endl;
			}
			
			//Realizamos la multiplicacion
			
			p_actual = 2;
			
			//Mientras el primo actual no supere el valor de K
			while(p_actual < k && !exito){
				p_potencia = p_actual;
				//Calculamos el valor de p_actual^e tal que no supere K
				while(p_potencia < k){
					p_potencia = p_potencia * p_actual;	//p_potencia = p_actual^{n+1}
				}
				//Dividimos para que p_potencia < k
				//mpz_divexact(p_potencia.get_mpz_t(), p_potencia.get_mpz_t(), p_actual.get_mpz_t());	
				//Multiplicamos dicho valor en el acumulador
				M = multiplicacionCurvaEliptica(Q, p_potencia, a, b, n, inv);
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
					}
				}
				//Pasamos al siguiente primo
				mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
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



mpz_class factorizacionCribaCuadratica(const mpz_class &n, const mpz_class &k, 
									   const unsigned int tam_tabla, bool debug){
    //Calculamos el numero de subtablas a crear
	const int tam_subtabla = 100000;
	const unsigned int num_tablas = tam_tabla / tam_subtabla + 1;
	vector<mpz_class> tabla_criba(tam_subtabla);
	vector<mpz_class> base_primos;
	vector<mpz_class> x_uniformes;
	vector<vector<bool>> factores_uniformes;
	vector<vector<int>> factores_uniformes_count;
	mpz_class p_actual = 3;
	mpz_class x_ini, entrada, r1, r2, sol1, sol2, it, p, q, aux;
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
			//if(debug) cout << p_actual << endl;
		}
		mpz_nextprime(p_actual.get_mpz_t(), p_actual.get_mpz_t());
	}
	if(debug)
		cout << "Tamanio base de primos: " << base_primos.size() << endl;
	
	
	
	for(int tabla_actual = 0; tabla_actual < num_tablas && x_uniformes.size() < base_primos.size()+50; tabla_actual++){
		if(debug) cout << "Tabla " << tabla_actual << endl;
		//Creamos la tabla de la criba y la tabla de factores
		mpz_sqrt(x_ini.get_mpz_t(), n.get_mpz_t());	//Guardamos en x_ini = sqrt(n)
		x_ini++;
		for(unsigned int i = 0; i < tam_subtabla; i++){
			entrada = i + x_ini + tabla_actual*tam_subtabla;
			//cout << "Entrada[" << entrada << "] = " << entrada << endl;
			entrada = entrada*entrada - n;
			tabla_criba[i] = entrada;
			//if(debug) cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << endl;
		}
		vector<vector<bool>> factores(tam_subtabla, vector<bool>(base_primos.size(), false));
		vector<vector<int>> factores_count(tam_subtabla, vector<int>(base_primos.size(), 0));
		//Calculamos las soluciones de f(x) = x^2 - n = 0 (mod p)
		//--> x^2 = n (mod p) --> x = sqrt(n) (mod p)
		//Una vez las tengamos, ejecutamos la criba
		for(unsigned int i = 0; i < base_primos.size(); i++){
			//Calculamos la raiz cuadrada de n
			//if(debug) cout << "Primo " << base_primos[i] << " " << endl;
			r1 = sqrtMod(n, base_primos[i]);
			
			//Realizamos la criba con la solucion positiva
			sol1 = (r1 - x_ini);
			mpz_mod(sol1.get_mpz_t(), sol1.get_mpz_t(), base_primos[i].get_mpz_t());
			//Calculamos, a partir de la primera posicion, la primera posicion en la tabla actual
			aux = (tam_subtabla*tabla_actual-sol1)/base_primos[i];
			sol1 = sol1 + base_primos[i]*aux - tam_subtabla*tabla_actual;
			if(sol1 < 0){
				sol1 += base_primos[i];
			}
			//if(debug) cout << "criba con " << sol1;
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
			//Calculamos, a partir de la primera posicion, la primera posicion en la tabla actual
			aux = (tam_subtabla*tabla_actual-sol2)/base_primos[i];
			sol2 = sol2 + base_primos[i]*aux - tam_subtabla*tabla_actual;
			if(sol2 < 0){
				sol2 += base_primos[i];
			}
			//Realizamos la criba con la solucion negativa
			if(sol1 != sol2){
				//if(debug) cout << " y criba con " << sol2;
				while(sol2 < tabla_criba.size()){
					tabla_criba[sol2.get_ui()] /= base_primos[i];
					factores[sol2.get_ui()][i] = !factores[sol2.get_ui()][i];
					factores_count[sol2.get_ui()][i] += 1;
					sol2 += base_primos[i];
				}
			}
			//if(debug) cout << endl;
		}
		
		for(unsigned int i = 0; i < tabla_criba.size(); i++){
			//Extraemos los restos uniformes (=1)
			if(tabla_criba[i] == 1){
				x_uniformes.push_back(i+x_ini+tabla_actual*tam_subtabla);
				factores_uniformes.push_back(factores[i]);
				factores_uniformes_count.push_back(factores_count[i]);
			}
			/*
			if(debug) {
				cout << "Tabla_criba["<<i<<"] = " << tabla_criba[i] << "\t[";
				
				for(unsigned int j = 0; j < base_primos.size(); j++){
					//cout << factores[i][j] << "\t";
					cout << factores_count[i][j] << "\t";
			}
			cout << "]" << endl;
			}
			*/
		}	
		if(debug) cout << "Restos uniformes encontrados: " << x_uniformes.size() << endl;
	}
	if(x_uniformes.size() > 0){
		if(debug){
			cout << "Restos k-uniformes: " << x_uniformes.size() << endl;
			cout << "Tamanio base de primos: "	<< base_primos.size() << endl;
		}
		/*
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
		*/
		//Trasponemos factores_uniformes y ejecutamos la eliminacion gaussiana
		factores_uniformes = transpose(factores_uniformes);
		
		factores_uniformes = eliminacionGaussiana(factores_uniformes);
		
		vector<vector<bool>> sols = encuentraSoluciones(factores_uniformes);
		if(debug) 
			cout << "Numero de soluciones: " << sols.size() << endl;
		vector<bool> sol;
		if(sols.size() > 1){
			cout << "Soluciones encontradas: " << sols.size() << endl;
			for(unsigned int n_sol = 0; n_sol < sols.size() && !exito; n_sol++){
				sol = sols[n_sol];
				if(debug) {
					cout << "Solucion \n[";
					for(unsigned int i = 0; i < sol.size(); i++){
						cout << sol[i] << "\t";
					}
					cout << "]" << endl;
				}
				//Si no es una solucion trivial
				if(any_of(sol.begin(), sol.end(), [](bool val) { return val; })){
					vector<int> factores_r(base_primos.size(), false);
					for(unsigned int i = 0; i < sol.size(); i++){
						if(sol[i]){
							for(unsigned int j = 0; j < base_primos.size(); j++){
								factores_r[j] += factores_uniformes_count[i][j];
							}
						}
					}
					
					//Calculamos la raiz cuadrada del numero dividiendo entre dos los exponentes
					for(unsigned int j = 0; j < base_primos.size(); j++){
						factores_r[j] /= 2;
					}
					
					//Calculamos la congruencia resultante x^2 = y^2 mod n
					mpz_class x_acum = 1;
					mpz_class r_acum = 1, r2_acum = 1;
					mpz_class x_actual;
					for(unsigned int i = 0; i < sol.size(); i++){
						if(sol[i] == 1){
							x_actual = x_uniformes[i];
							/*
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
							*/
							x_acum *= x_actual;
						}
					}
					
					for(unsigned int j = 0; j < factores_r.size(); j++){
						for(int k = 0; k < factores_r[j]; k++)
							r_acum *= base_primos[j];
					}
					
					if(debug) cout << "Congruencia " << n_sol << endl;
					
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
				else{
					if(debug)
						cout << "Saltando solucion trivial..." << endl;
				}
			}
		}
		else{
			if(debug) cout << "El sistema no tiene solucion no trivial :(" << endl;
			p = n;
		}
	}
	else{
		if(debug) cout << "No se ha encontrado ningun resto k-uniforme :(" << endl;
		p = n;
	}
	return p;
}

bool ataqueWiener(mpz_class &d, mpz_class &p, mpz_class &q, 
				  const mpz_class e, const mpz_class n, bool debug){
	
	mpz_class convergente_num, convergente_den;
	mpz_class k1 = 1, k2 = 0, d1 = 0, d2 = 1, k;
	mpz_class phi, r, a, b, c;
	mpz_class aux;
	bool exito = false, existe_solucion = false;
	vector<mpz_class> cf = cocientesFraccionContinua(e, n);
	if(debug) cout << "e = " << e << endl;
	if(debug) cout << "n = " << n << endl;
	for(unsigned int i = 0; i < cf.size() && !exito; i++){
		//Calculamos la reducida i-esima
		d = cf[i] * d1 + d2;
        k = cf[i] * k1 + k2;
		if(debug) cout << "Reducida " << i << " = " << k << "/" << d << endl;
        d2 = d1;
        d1 = d;
        k2 = k1;
	    k1 = k;
	    
		//Si divide entre 0, pasamos al siguiente convergente
		if(k != 0){
			//Calculamos (e*d-1)/k 
			aux = e*d-1;
			mpz_fdiv_qr(phi.get_mpz_t(), r.get_mpz_t(), aux.get_mpz_t(), k.get_mpz_t());
			if(debug){
				cout << "resto = " << r << endl;
				cout << "Potencial phi = " << phi << endl;
			}
			//Si (e*d-1)/k no es una division exacta, pasamos a la siguiente convergente
			if(r == 0){
				//Resolvemos la ecuacion x^2 - (N - phi + 1) x + N = 0 --> Raices = p,q
				a = 1;
				b = -(n - phi + 1);
				c = n;
				existe_solucion = resuelveEcuacionCuadratica(p, q, a, b, c);
				if(debug){
					cout << "existe solucion = " << existe_solucion << endl;
					cout << "p*q = " << p*q << endl;
				}
				if (existe_solucion && p * q == n) {
				    exito = true;
				}
			}
		}
	}
	
	if(!exito){
		d = 1;
		p = n;
		q = 1;
	}
	
	return exito;
}

