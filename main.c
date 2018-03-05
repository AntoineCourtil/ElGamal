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



/**
* Fonction qui effectue nbIteration de test d'Eucide
**/

void testEuclide(int nbIteration, int affichage){
	mpz_t g, p, u, v;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("euclideResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}

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

		//Affichage resultats
		if(affichage){
			printf("\n\n\nResultat : \n");
			printauplusbv(g, u, p, v, pgcd);
		}



		//Ecriture dans le fichier results

		char* gString = mpz_get_str(NULL, 10, g);
		char* uString = mpz_get_str(NULL, 10, u);
		char* pString = mpz_get_str(NULL, 10, p);
		char* vString = mpz_get_str(NULL, 10, v);
		char* pgcdString = mpz_get_str(NULL, 10, pgcd);

		fprintf(file, "%s\n     *(\n%s\n     )\n     +\n%s\n     *(%s\n    )\n     =\n%s\n\n\n", gString, uString, pString, vString, pgcdString);

	}

	fclose(file);
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


/**
* Fonction qui effectue nbIteration de test de expMod
**/

void testExpMod(int nbIteration, int affichage){
	mpz_t g, p, a;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("expModResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}

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

		//Affichage resultats
		if(affichage){
			printf("\nResultats g=%d : res = ", i);
			mpz_out_str(NULL, 10, res);
			printf(" et res2 = ");
			mpz_out_str(NULL, 10, res2);
		}



		//Ecriture dans le fichier results

		char* resString = mpz_get_str(NULL, 10, res);
		char* res2String = mpz_get_str(NULL, 10, res2);

		fprintf(file, "     g : %d\nres : %s\nverifRes : %s\n\n\n", i, resString, res2String);


		//Liberation de mémoire
		mpz_clears(res, res2, (void *) NULL);


	}

	fclose(file);
}

/**
* Calcul pour x un grand nombre random en fonction de la borne p
**/

void getRandomMpzt(mpz_t x, mpz_t p){


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
}



/**
* Génére la clé publique publicKey
* en fonction de p et g
**/

void keygen(mpz_t p, mpz_t g, mpz_t publicKey, int affichage){

	
	//Init de x
	mpz_t x;
	mpz_init(x);

	//Calcul du random pour x
	getRandomMpzt(x, p);


	//Calcul de la Public Key
	expMod(p, g, x, publicKey);	


	//Affichage des résultats
	if(affichage == 1){
		printPublicKey(p, g, x, publicKey);

		//Verification
		mpz_powm(publicKey, g, x, p);

		printf("\n\n          Verification :\n");
		printPublicKey(p, g, x, publicKey);

		printf("\n\n___________________________________\n\n");
	}

}


/**
* Fonction qui effectue nbIteration de test de expMod
**/

void testKeygen(int nbIteration){

	mpz_t p, g, publicKey;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("keygenResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}

	for(i=1;i<=nbIteration;i++){

		mpz_init(p);
		mpz_init(g);
		mpz_init(publicKey);

		//Initialisation de g
		mpz_set_ui(g, i);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);

		//Calcul de la publicKey
		keygen(p, g, publicKey, 0);


		//Ecriture dans le fichier results

		char* pString = mpz_get_str(NULL, 10, p);
		char* gString = mpz_get_str(NULL, 10, g);
		char* publicKeyString = mpz_get_str(NULL, 10, publicKey);

		fprintf(file, "     p : %s\n     g : %s\npublicKey : %s\n\n\n", pString, gString, publicKeyString);


		//Liberation de mémoire
		mpz_clears(p, g, publicKey, (void *) NULL);
	}

	fclose(file);
}


/**
* Chiffre un message m en fonction de la clé publique publicKey
**/

void encrypt(mpz_t publicKey, mpz_t p, mpz_t g, mpz_t m, mpz_t C, mpz_t B, mpz_t x, int recalculRandom){
	

	if(recalculRandom){
		//Calcul du random pour x
		getRandomMpzt(x, p);
	}

	//Init de y mod p
	mpz_t y;
	mpz_init(y);

	//Calcul de y
	expMod(p, publicKey, x, y);

	//y = y mod p
	mpz_mod(y, y, p);

	//C = m * y mod p
	mpz_mul(C, m, y);

	//C = C mod p
	mpz_mod(C, C, p);

	//B = g^x mod p
	expMod(p, g, x, B);

}


/**
* Fonction qui effectue nbIteration de test de encrypt
**/

void testEncrypt(int nbIteration, int affichage){
	mpz_t p, g, publicKey, C, B, m, x;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("encryptResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}

	for(i=1;i<=nbIteration;i++){

		mpz_init(p);
		mpz_init(g);
		mpz_init(publicKey);
		mpz_init(C);
		mpz_init(B);
		mpz_init(m);

		//Initialisation de g
		mpz_set_ui(g, 2);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);

		//Calcul de la publicKey
		keygen(p, g, publicKey, 0);

		//Initialisation de m
		mpz_set_ui(m, i);

		//Chiffrement de m
		encrypt(publicKey, p, g, m, C, B, x, 1);


		//Affichage des résultats
		if(affichage){
			printf("\nm : %d", i);
			printf("\nC : ");
			mpz_out_str(NULL, 10, C);
			printf("\nB : ");
			mpz_out_str(NULL, 10, B);
			printf("\n\n");
		}


		//Ecriture dans le fichier results

		char* CString = mpz_get_str(NULL, 10, C);
		char* BString = mpz_get_str(NULL, 10, B);

		fprintf(file, "m : %d\nC : %s\nB : %s\n\n\n", i, CString, BString);


		//Liberation de mémoire
		mpz_clears(p, g, publicKey, C, B, m, (void *) NULL);
	}

	fclose(file);
}


/**
* Déchiffre un message m en fonction de la clé secret secretKey
**/

void decrypt(mpz_t C, mpz_t B, mpz_t secretKey, mpz_t m, mpz_t p){

	//Init de D pour le déchriffrement
	mpz_t D;
	mpz_init(D);

	//Init de u et v pour euclide
	mpz_t u, v;
	mpz_init(u);
	mpz_init(v);

	//Init de pgcd pour le déchriffrement
	mpz_t pgcd;
	mpz_init(pgcd);


	//D = B^secretKey mod p
	expMod(p, B, secretKey, D);


	//pgcd de D et p
	//avec u = D^-1
	euclide(D, u, p, v, pgcd);


	//m = C × D^−1 mod p.
	mpz_mul(m, C, u);


	//m = m mod p
	mpz_mod(m, m, p);


	//Liberation de mémoire
	mpz_clears(D, pgcd, u, v, (void *) NULL);
}


/**
* Fonction qui effectue nbIteration de test de decrypt
**/

void testDecrypt(int nbIteration, int affichage){

	mpz_t p, g, publicKey, C, B, m, secretKey, mDecrypt;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("decryptResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}

	for(i=1;i<=nbIteration;i++){

		mpz_init(p);
		mpz_init(g);
		mpz_init(publicKey);
		mpz_init(C);
		mpz_init(B);
		mpz_init(m);
		mpz_init(secretKey);
		mpz_init(mDecrypt);

		//Initialisation de g
		mpz_set_ui(g, 2);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);

		//Calcul de la publicKey
		keygen(p, g, publicKey, 0);

		//Initialisation de m
		mpz_set_ui(m, i);

		//Chiffrement de m
		encrypt(publicKey, p, g, m, C, B, secretKey, 1);

		//Dechiffrement de m
		decrypt(C, B, secretKey, mDecrypt, p);


		//Affichage des résultats
		if(affichage){
			printf("\nm : %d", i);
			printf("\n     C : ");
			mpz_out_str(NULL, 10, C);
			printf("\n     B : ");
			mpz_out_str(NULL, 10, B);
			printf("\n     mDecrypt : ");
			mpz_out_str(NULL, 10, mDecrypt);
			printf("\n\n");
		}



		//Ecriture dans le fichier results

		char* CString = mpz_get_str(NULL, 10, C);
		char* BString = mpz_get_str(NULL, 10, B);
		char* mDecryptString = mpz_get_str(NULL, 10, mDecrypt);

		fprintf(file, "m : %d\nC : %s\nB : %s\nmDecrypt : %s\n\n\n", i, CString, BString, mDecryptString);



		//Liberation de mémoire
		mpz_clears(p, g, publicKey, C, B, m, mDecrypt, secretKey, (void *) NULL);
	}


	fclose(file);
}



/**
* Vérifie la propriété homomorphique de ElGamal 
**/

void homomorphie(mpz_t p, mpz_t g, mpz_t publicKey, mpz_t secretKey, int affichage, FILE* file){

	//Init des variables à calculer
	mpz_t C, C1, C2, B, B1, B2;

	mpz_init(C);
	mpz_init(C1);
	mpz_init(C1);
	mpz_init(B);
	mpz_init(B1);
	mpz_init(B2);



	//Init des messages m
	mpz_t m, m1, m2, m1xm2;

	mpz_init(m);
	mpz_init(m1);
	mpz_init(m2);
	mpz_init(m1xm2);


	//Calcul des messages m1 et m2
	mpz_set_ui(m1, rand());
	mpz_set_ui(m2, rand());
	mpz_mul(m1xm2, m1, m2);


	//Calcul de C1 et B1
	encrypt(publicKey, p, g, m1, C1, B1, secretKey, 1);

	//Calcul de C2 et B2
	encrypt(publicKey, p, g, m2, C2, B2, secretKey, 0);


	//C = C1 * C2
	mpz_mul(C, C1, C2);

	//B = B1 * B2
	mpz_mul(B, B1, B2);


	//Dechiffrement de C et B dans m
	decrypt(C, B, secretKey, m, p);

	//m = m mod p
	mpz_mod(m, m, p);


	//Affichage des résultats
	if(affichage){
		printf("\n     m1 : ");
		mpz_out_str(NULL, 10, m1);
		printf("\n     m2 : ");
		mpz_out_str(NULL, 10, m2);
		printf("\n     m1xm2 : ");
		mpz_out_str(NULL, 10, m1xm2);
		printf("\n     m    : ");
		mpz_out_str(NULL, 10, m);
		printf("\n\n");
	}

	//Ecriture dans le fichier result

	char* m1String = mpz_get_str(NULL, 10, m1);
	char* m2String = mpz_get_str(NULL, 10, m2);
	char* m1xm2String = mpz_get_str(NULL, 10, m1xm2);
	char* mString = mpz_get_str(NULL, 10, m);

	fprintf(file, "m1 : %s\nm2 : %s\nm1xm2 : %s\nm decrypt : %s\n\n\n", m1String, m2String, m1xm2String, mString);

	
}


/**
* Fonction qui effectue nbIteration de test de decrypt
**/

void testHomomorphie(int nbIteration){

	mpz_t p, g, publicKey, secretKey;

	int i;

	//Init fichier de sortie
	FILE* file = fopen("homomorphieResults.txt", "w");

	if(file == NULL){
		printf("Erreur d'ouverture de fichier\n");
	}


	

	for(i=1;i<=nbIteration;i++){

		mpz_init(p);
		mpz_init(g);
		mpz_init(publicKey);
		mpz_init(secretKey);

		//Initialisation de g
		mpz_set_ui(g, 2);

		//Initialisation de p
		mpz_init_set_str(p, P_HEXVALUE, 16);

		//Calcul de la publicKey
		keygen(p, g, publicKey, 0);

		//Test de la propriété homomorphique
		homomorphie(p, g, publicKey, secretKey, 0, file);

		//Liberation de mémoire
		mpz_clears(p, g, publicKey, secretKey, (void *) NULL);
	}



	fclose(file);
}


 
int main()
{
	testEuclide(10000, 0);

	testExpMod(100, 0);

	testKeygen(100);

	testEncrypt(100, 0);

	testDecrypt(100, 0);

	testHomomorphie(100);

	return 0; 
}