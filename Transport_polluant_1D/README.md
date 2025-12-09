# Transport d'un polluant 1D

Ce projet contient des scripts Python permettant d’étudier la propagation d’un polluant dans un fluide en écoulement uniforme selon l’axe \(x\), sans diffusion.


## Description

Le projet est organisé en un script principal :

**transport.py** : définition des fonctions pour résoudre l’équation de transport :

- `thomas(A,D)` : résolution d’un système tridiagonal par l’algorithme de Thomas.
- `Courant(L,V,beta,t,ntime,N)` : schéma de Courant pour l’équation de transport.
- `LaxWendroff(L,V,beta,t,ntime,N)` : schéma de Lax-Wendroff.
- `LeapFrog(L,V,beta,t,ntime,N)` : schéma Leap-Frog.
- `sol_ex(x,t,V,L,beta)` : solution analytique de l’équation de transport.

Le script principal compare les solutions numériques aux solutions analytiques pour différentes méthodes et valeurs du paramètre $r = V \Delta t / \Delta x$.



## Fonctionnalités

- Calcul de la concentration \(C(x,t)\) à différents instants.
- Comparaison de la solution numérique avec la solution analytique.
- Étude de l’effet de différents schémas numériques :
  - Courant
  - Lax-Wendroff
  - Leap-Frog
- Visualisation graphique des résultats avec sauvegarde des figures.



## Outils utilisés

- Python
- Numpy
- Matplotlib



## Résultats

- Graphiques des solutions numériques et exactes pour différents schémas et temps.
- Étude de la diffusion et dispersion numérique selon le schéma choisi et le paramètre \(r\).



## Auteur

Mouhalhal Moussaab




