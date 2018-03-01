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
      printf("a * u   +   p * v   = pgcd(a,p)\n");

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

		printf("u : ");
		mpz_out_str(NULL, 10, u);
		printf("   -   v : ");
		mpz_out_str(NULL, 10, v);
		printf("\n");




		//Decalage de valeurs
		mpz_set(res0, res1);
		mpz_set(res1, r);

		mpz_set(u0, u1);
		mpz_set(u1, u);

		mpz_set(v0, v1);
		mpz_set(v1, v);
	}
}
 
int main()
{
	mpz_t a, p, u, v, g;

	mpz_init(a);
	mpz_init(p);
	mpz_init(u);
	mpz_init(v);
	mpz_init(g);


	//Init pgcd
	mpz_t pgcd;
	mpz_init(pgcd);

	//Initialisation de g
	mpz_set_ui(g, 2);

	//Initialisation de p
	//mpz_init_set_str(p, P_HEXVALUE, 16);


	mpz_set_ui(a, 1234);
	mpz_set_ui(p, 5678);



	euclide(a, u, p, v, pgcd);

	printf("\n\n\nResultat : \n");
	printauplusbv(a, u, p, v, pgcd);

    return 0; 
}