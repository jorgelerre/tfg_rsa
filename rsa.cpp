#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>

using namespace std;

bool strongPseudoprime(const mpz_t p, const mpz_t base, bool debug = true){
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
bool BigInt::strongPseudoprime(const BigInt& base, bool debug) const{
	bool is_pseudoprime = false;
	BigInt mcd, mcm, u0, v0;
	BigInt q, r;
	BigInt b;
	BigInt zero = "0x0", one = "0x1", two = "0x2";
	//Comprobamos que sean coprimos y que el primo que comprobamos sea mayor que 1
	this->EEA(base, mcd, mcm, u0, v0);
	if(debug)
		cout << "mcd de a y la base = " << mcd << endl;
	if(mcd == one && one < (*this)){
		//Calculamos s,t tales que n-1 = 2^s*t con t impar
		BigInt t = (*this) - one, s = zero;
		t.scholarDivision(two, q, r);
		while (r == zero){
			t = q;
			s = s + one;
			t.scholarDivision(BigInt(two), q, r);
		}
		if(debug){
			cout << "s = " << s << endl;
			cout << "t = " << t << endl;
			cout << "2^s*t = " << two.quickModExp(s,*this)*t << endl;	//el modulo es mayor que 2^s
		}
		
		b = base.quickModExp(t, *this);
		
		if(b == one || b == (*this) - one){
			is_pseudoprime = true;
		}
		for(BigInt r = one; r < s - one && !is_pseudoprime; r = r + one){
			b = (b*b) % (*this);
			if(b == (*this) - one){
				is_pseudoprime = true;
			}
		}
	}
	return is_pseudoprime;
}

bool BigInt::millerRabinTest(int k, bool debug){
	bool no_fail = true;
	BigInt rn;
	for(int i = 0; i < k && no_fail; i++){
		rn = randomBigInt(*this);
		if(debug){
			cout << "-------------i = " << i << "-------------" << endl;
			cout << "Numero aleatorio generado = " << rn << endl;
		}
		if(strongPseudoprime(rn, debug) == false){
			no_fail = false;
		}
	}
	return no_fail;
}
*/

/*
 * Genera un numero aleatorio 
 * @p Numero entero a analizar
 * @it Numero de iteraciones que aplicar la prueba
*/
void generate_random_limb(mp_limb_t &random_limb, gmp_randstate_t state) {
    mpz_t random_num;
    mpz_init(random_num);

    // Genera un número aleatorio de tamaño máximo de un mp_limb_t
    mpz_urandomb(random_num, state, GMP_LIMB_BITS);

    // Convierte el número generado a mp_limb_t
    random_limb = mpz_getlimbn(random_num, 0);

    // Limpieza
    mpz_clear(random_num);
}

/*
 * Aplica el test de Miller Rabin sobre un primo cierto numero de veces
 * @p Numero entero a analizar
 * @it Numero de iteraciones que aplicar la prueba
*/
bool millerRabin(mpz_t p, int it){
	bool no_fail = true;
	mp_limb_t a;
	for(int i = 0; i < it && no_fail; i++){
		
	}
	
	return no_fail;
}


// Función para generar un número primo de 'bits' bits
void generate_prime(mpz_t prime, int bits) {
	gmp_randstate_t gmp_randstate;
	
    mpz_rrandomb(prime, gmp_randstate, bits); // Genera un número aleatorio de 'bits' bits
    mpz_nextprime(prime, prime); // Encuentra el siguiente número primo
}

// Función principal para generar claves RSA
void generate_rsa_key(int bits, mpz_t n, mpz_t e, mpz_t d) {

}

int main() {
	bool pseudoprimo;
	mpz_t p, base;
	mpz_inits(p, base, nullptr);
	mpz_set_ui(p, 501);
	mpz_set_ui(base, 2);
	pseudoprimo = strongPseudoprime(p, base, true);
	cout << pseudoprimo << endl;
	mpz_clears(p, base, nullptr);
    return 0;
}

