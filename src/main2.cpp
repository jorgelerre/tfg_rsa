#include "rsa.h"
#include "rsa_attacks.h"
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

void muestraMenuPrincipal() {
    cout << "========== Menu ===========\n";
    cout << "1. Crear clave RSA\n";
    cout << "2. Cifrar mensaje\n";
    cout << "3. Descifrar mensaje\n";
    cout << "4. Atacar clave\n";
    cout << "5. Salir\n";
    cout << "===========================\n";
    cout << "Seleccione una opción: ";
}

void muestraMenu1() {
    cout << "========= Menu 1 ==========\n";
    cout << "1. Primos aleatorios\n";
    cout << "2. Primos robustos\n";
    cout << "3. Volver\n";
    cout << "===========================\n";
    cout << "Seleccione una opción: ";
}


void muestraMenu4() {
    cout << "========= Menu 2 ==========\n";
    cout << "1. Ataque de Fermat\n";
    cout << "2. Ataque de Kraitchik\n";
    cout << "3. Ataque rho de Pollard\n";
    cout << "4. Ataque p-1\n";
    cout << "5. Ataque de factorizacion por curvas elipticas\n";
    cout << "6. Ataque de factorizacion por criba cuadratica\n";
    cout << "7. Ataque de Wiener\n";
    cout << "8. Volver\n";
    cout << "===========================\n";
    cout << "Seleccione una opción: ";
}

void muestraMenu3() {
	
	
}



int main() {
	bool seguir = true, debug = true, clave_iniciada = false;
	mpz_class n, e, d, p, q;
    int opcion, opcion2;
    unsigned int tam_n, tam_tabla, att;
    string input;
    mpz_class mensaje;
    mpz_class k;
    
	//Inicializacion de numeros aleatorios
	gmp_randstate_t state;
	gmp_randinit_default(state);
    gmp_randseed_ui(state, static_cast<unsigned long>(time(nullptr)));
    
    
    //Programa principal - Menu
	while (seguir) {
        muestraMenuPrincipal();
        cin >> opcion;
        
        switch (opcion) {
		    case 1:
		    	muestraMenu1();
		    	cin >> opcion2;
		        while(opcion2 < 1 || opcion2 > 3){
			        cout << "Opción inválida. Por favor, intentalo de nuevo.\n" << endl;
		        	cin >> opcion2;
		        }
		        //Si la opcion no es 3, pedimos el tamanio deseado de clave
		        if(opcion2 != 3){
		        	cout << "Introduce el numero de bits deseado para n: " ;
		        	cin >> tam_n;
		        	while(tam_n < 0){
		        		cout << "El numero debe ser mayor que 5. Intentalo de nuevo: " << endl;
		        		cin >> tam_n;
		        	}
		        	if(opcion2 == 1){
		        		generate_rsa_key(n, e, d, tam_n, state, false, debug);
		        	}
		        	else{
		        		generate_rsa_key(n, e, d, tam_n, state, true, debug);
		        	}
		        	cout << "----Clave generada----" << endl;
		        	cout << "n = " << n << endl;
		        	cout << "e = " << e << endl;
		        	cout << "d = " << d << endl;
		        	cout << "tamanio n = " << mpz_sizeinbase(n.get_mpz_t(), 2) << endl;
		        	clave_iniciada = true;
				}
		        break;
		    case 2:
		        if(clave_iniciada){
		        	cout << "Introduce el numero a cifrar: ";
		        	cin >> mensaje;
		        	cout << "Mensaje = " << mensaje << endl;
		        	cout << "Cifrado = " << cifra_RSA(mensaje, e, n) << endl;
		        }
		        else{
		        	cout << "Necesitas inicializar la clave para usar esta opcion.\n";
		        }
		        break;
		    case 3:
		        if(clave_iniciada){
		        	cout << "Introduce el numero a descifrar: ";
		        	cin >> mensaje;
		        	cout << "Cifrado = " << mensaje << endl;
		        	cout << "Mensaje original = " << descifra_RSA(mensaje, d, n) << endl;
		        }
		        else{
		        	cout << "Necesitas inicializar la clave para usar esta opcion.\n";
		        }
		        break;
		    case 4:
		        if(clave_iniciada){
		        	muestraMenu4();
		        	cin >> opcion2;
		        	while(opcion2 < 1 || opcion2 > 8){
					    cout << "Opción inválida. Por favor, intentalo de nuevo.\n" << endl;
				    	cin >> opcion2;
				    }
				    switch(opcion2){
				    	case 1:	//Factorizacion de Fermat
				    		cout << "Factorizacion de Fermat\n";
				    		p = ataqueFermat(n, state);
				    		break;
				    	case 2: //Factorizacion de Kraitchik
				    		cout << "Factorizacion de Kraitchik\n";
				    		p = ataqueKraitchik(n, state);
				    		break;
				    	case 3:	//Rho de Pollard
				        	cout << "Rho de Pollard\n";
				        	p = rhoPollard(n, state);
				        	break;
				        case 4:	//p-1 de Pollard
							cout << "Ataque p-1\n";
							break;
						case 5: //Ataque de factorizacion por curvas elipticas
							cout << "Ataque de factorizacion por curvas elipticas\n";
							cout << "Introduce un valor de k: ";
							cin >> k;
							cout << "k = " << k << endl;
							cout << "Introduce el numero de intentos a realizar: ";
							cin >> att;
							cout << "att = " << k << endl;
							p = factorizacionCurvasElipticas(n, state, k, att);
							break;
						case 6:	//Ataque de factorizacion por criba cuadratica
							cout << "Ataque de factorizacion por criba cuadratica\n";
							cout << "Introduce un valor de k: ";
							cin >> k;
							cout << "k = " << k << endl;
							cout << "Introduce un valor de tam_tabla: ";
							cin >> tam_tabla;
							cout << "tam_tabla = " << tam_tabla << endl;
							p = factorizacionCribaCuadratica(n, k, tam_tabla, true);
							break;
						case 7: //Ataque de Wiener
							cout << "Ataque de Wiener\n";
							break;
						case 8: //Volver
							break;
				    }
				    if(opcion2 != 8){
						q = n / p;
						cout << "n = " << n << endl;
						cout << "p = " << p << endl;
						cout << "q = " << q << endl;
					}
		        }
		        else{
		        	cout << "Necesitas inicializar la clave para usar esta opcion.\n";
		        }
		        break;
		    case 5:
		        cout << "Saliendo del programa...\n";
		        seguir = false;
		        break;
		    default:
		        cout << "Opción inválida. Por favor, intenta de nuevo.\n";
		        break;
    	}
    }
	
	
    return 0;
}

