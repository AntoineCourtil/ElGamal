# Rapport 

## Question 1. Quel langage de programmation avez-vous choisi ?

Comme le conseillait le sujet, j'ai choisi le langage C dans le but d'intégrer la bibliothèque GMP afin de pouvoir effectuer les opérations nécessaires sur les grands nombres.

___

## Question 2. En vous aidant d’internet, donnez la définition d’un nombre aléatoire cryptographiquement sûr.

Un générateur de nombres aléatoires, est un dispositif capable de produire une séquence de nombres pour lesquels, il n'existe aucun lien déterministe entre un nombre et ses prédécesseurs, de façon que cette séquence puisse être appelée « suite de nombres aléatoires ».

On évalue la qualité d'un générateur de nombre aléatoire en fonction de sa période et de sa discrépance.

- Sa période est le nombre d'itérations à partir duquel le générateur recalcule la même suite de nombre aléatoires. Pour être considéré comme convenable aujourd’hui, un générateur se doit de posséder une période d’au moins 2^100.

- Sa discrépance mesure l'écart maximal de la suite donnée avec l'équirépartition.


___

## Question 3. Implémentez la fonction Euclide

___

## Question 4. Implémentez la fonction ExpMod

J'ai choisi de créer les variables :
- newA (nouvelle valeur de a en fonction de g pair ou impair)
- g2 (g^2) 
Afin de ne pas écraser les valeur de a et g initiale.


___

## Question 5. Implémentez ces trois fonctions

- Keygen :

- Encrypt :
	Etant donné que le state de mon random est lié au time de la machine, et qu'on exécute plusieurs fois en même temps, le B est le même.

- Decrypt :


___

## Question 6. Vérifiez la propriété homomorphique en chiffrant deux messages différents


___

### Sources :

http://math.univ-lyon1.fr/~jberard/genunif-www.pdf