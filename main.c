#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

#define P_HEXVALUE "FFFFFFFF FFFFFFFF C90FDAA2 2168C234 C4C6628B 80DC1CD1 29024E08 8A67CC74 020BBEA6 3B139B22 514A0879 8E3404DD EF9519B3 CD3A431B 302B0A6D F25F1437 4FE1356D 6D51C245 E485B576 625E7EC6 F44C42E9 A637ED6B 0BFF5CB6 F406B7ED EE386BFB 5A899FA5 AE9F2411 7C4B1FE6 49286651 ECE65381 FFFFFFFF FFFFFFFF"



/**
* Affiche au + bv = p
**/

void printauplusbv(mpz_t a, mpz_t u, mpz_t b, mpz_t v, mpz_t p) 
{
      //printf("a * u   +   p * v   = pgcd(a,p)\n");

      mpz_out_str(NULL, 10, a);
      printf(" * (");
      mpz_out_str(NULL, 10, u);
      printf(") + ");
      mpz_out_str(NULL, 10, b);
      printf(" * (");
      mpz_out_str(NULL, 10, v);
      printf(") = ");
      mpz_out_str(NULL, 10, p);
      printf("\n");
}


void printPublicKey(mpz_t p, mpz_t g, mpz_t x, mpz_t publicKey){
	printf("\nPublicKey :\n     p = ");
	mpz_out_str(NULL, 10, p);
	printf("\n     g = ");
	mpz_out_str(NULL, 10, g);
	printf("\n     x = ");
	mpz_out_str(NULL, 10, x);
	printf("\n     publicKey = ");
	mpz_out_str(NULL, 10, publicKey);
	printf("\n");
}





/**
*
* Calcul les coefficients de Bézout tel que
* a.u + p.v = pgcd(a,p)
*
**/

void euclide(mpz_t a, mpz_t u0, mpz_t p, mpz_t v0, mpz_t res0)
{

	//Copie de a et p pour calculs
	mpz_t r, res1;
	mpz_init(res1);
	mpz_init(r);

	mpz_set(res0, a);
	mpz_set(res1, p);

	//Init des premiers u
	mpz_t u, u1;
	mpz_init(u0);
	mpz_init(u1);
	mpz_init(u);

    mpz_set_ui(u0, 1);
    mpz_set_ui(u1, 0);

    //Init des premiers v
	mpz_t v, v1;
	mpz_init(v0);
	mpz_init(v1);
	mpz_init(v);

    mpz_set_ui(v0, 0);
    mpz_set_ui(v1, 1);

    //Init quotient
    mpz_t q;
    mpz_init(q);




    //Tant que le signe de res1 est positif
    while(mpz_sgn(res1) == 1) {


    	//printauplusbv(a, u, p, v, res0);

		//Var q = res0/res1   &&   Var r = res0%res1
		mpz_tdiv_qr(q, r, res0, res1);

		//Init des vars u, v
		mpz_set(u, u0);
		mpz_set(v, v0);


		//Calcul des coefficients stockés dans u et v
		mpz_submul(u, q, u1); // u = u-q*u1
		mpz_submul(v, q, v1); // v = v-q*v1

		/*printf("u : ");
		mpz_out_str(NULL, 10, u);
		printf("   -   v : ");
		mpz_out_str(NULL, 10, v);
		printf("\n");*/


		//Decalage de valeurs
		mpz_set(res0, res1);
		mpz_set(res1, r);

		mpz_set(u0, u1);
		mpz_set(u1, u);

		mpz_set(v0, v1);
		mpz_set(v1, v);
	}
}



void testEuclide(int nbIteration){
	mpz_t g, p, u, v;

	int i;

	for(i=1; i<=nbIteration; i++){

		mpz_init(g);
		mpz_init(p);
		mpz_init(u);
		mpz_init(v);


		//Init pgcd
		mpz_t pgcd;
		mpz_init(pgcd);

		//Initialisation de g
		mpz_set_ui(g, i);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);


		euclide(g, u, p, v, pgcd);

		printf("\n\n\nResultat : \n");
		printauplusbv(g, u, p, v, pgcd);
	}
}


/**
* Calcule l’exponentiation binaire de
* r = g^a mod p
**/

void expMod(mpz_t p, mpz_t g, mpz_t a, mpz_t res){

	//Reste de a/2 afin de savoir s'il est pair ou non
	mpz_t aDivBy2;
	mpz_init(aDivBy2);

	//Resultat de la nouvelle valeur de a pour continuer les calculs
	mpz_t newA;
	mpz_init(newA);

	//Resultat de g^2
	mpz_t g2;
	mpz_init(g2);

	//printf("\n\na : ");
	//mpz_out_str(NULL, 10, a);

	//si a=1
	if(mpz_cmp_ui(a, 1) == 0){
		mpz_mod(res, g, p);
	}
	else{

		//Calcul du reste de a/2
		mpz_mod_ui(aDivBy2, a, 2);

		//Si a est pair
		if(mpz_cmp_ui(aDivBy2, 0) == 0){

			//Calcul de a/2
			mpz_tdiv_q_ui(newA, a, 2);
			mpz_mod(newA, newA, p);

			//Calcul de g^2
			mpz_pow_ui(g2, g, 2);
			mpz_mod(g2, g2, p);

			//Calcul de res mod p
			mpz_mod(res, res, p);

			//Calcul de la prochaine recursion
			expMod(p, g2, newA, res);
		}

		//Si a est impair
		else{

			//Si a>2 est impair
			if(mpz_cmp_ui(a, 2) > 0){

				//Calcul de (a-1)/2
				mpz_sub_ui(newA, a, 1);
				mpz_tdiv_q_ui(newA, newA, 2);
				mpz_mod(newA, newA, p);

				//Calcul de g^2
				mpz_pow_ui(g2, g, 2);
				mpz_mod(g2, g2, p);

				//Calcul de la prochaine recursion
				expMod(p, g2, newA, res);

				//Calcul de res mod p
				mpz_mod(res, res, p);

				//r = r*g
				mpz_mul(res, res, g);

				//Calcul de res mod p
				mpz_mod(res, res, p);				
			}
		}
	}

	//Liberation de mémoire
	mpz_clears(aDivBy2, g2, newA, (void *) NULL);
}


void testExpMod(int nbIteration){
	mpz_t g, p, a;

	int i;

	for(i=1; i<=nbIteration; i++){

		mpz_init(g);
		mpz_init(p);
		mpz_init(a);


		//Init resultat res
		mpz_t res;
		mpz_init(res);

		//Init resultat qui doit etre trouvé res2
		mpz_t res2;
		mpz_init(res2);

		//Initialisation de g
		mpz_set_ui(g, i);

		//Initialisation de a
		mpz_set_ui(a, 2);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);


		expMod(p, g, a, res);
		mpz_powm(res2, g, a, p);

		printf("\nResultats g=%d : res = ", i);
		mpz_out_str(NULL, 10, res);
		printf(" et res2 = ");
		mpz_out_str(NULL, 10, res2);


		//Liberation de mémoire
		mpz_clears(res, res2, (void *) NULL);


	}
}



void keygen(mpz_t p, mpz_t g, mpz_t publicKey){

	//Init de x
	mpz_t x;
	mpz_init(x);


	//Init de la borne max pour le random de x 
	mpz_t borne;
	mpz_init(borne);

	//borne = p-1
	mpz_sub_ui(borne, p, 1);


	//Init seed pour le randseed
	long seed;

	//Ici le seed sera lié au time
	time(&seed);

	//Init state pour créer le random
	gmp_randstate_t state;
	gmp_randinit_default(state);

	//Création du state
	gmp_randseed_ui(state, seed);

	//Calcul du random pour x
	mpz_urandomm(x, state, borne);


	//Calcul de la Public Key
	expMod(p, g, x, publicKey);

	//Affichage des résultats
	printPublicKey(p, g, x, publicKey);


	//Verification
	mpz_powm(publicKey, g, x, p);

	//Affichage des résultats
	printf("\nVerification :");
	printPublicKey(p, g, x, publicKey);

}


void testKeygen(){

	mpz_t p, g, publicKey;

	mpz_init(p);
	mpz_init(g);
	mpz_init(publicKey);

	//Initialisation de g
	mpz_set_ui(g, 2);

	//Initialisation de p
	mpz_init_set_str(p, P_HEXVALUE, 16);

	keygen(p, g, publicKey);


}


 
int main()
{
	//testEuclide(10000);

	testExpMod(100);

	testKeygen();

    return 0; 
}