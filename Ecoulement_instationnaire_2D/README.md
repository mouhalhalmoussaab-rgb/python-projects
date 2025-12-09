# Écoulement instationnaire 2D

Ce projet contient des scripts Python permettant d’étudier l’évolution de la vitesse d’un écoulement dans une conduite rectangulaire à partir de l’équation de diffusion instationnaire en 2D.



## Description

Le projet est organisé en un script principal :

**code.py** : définition des fonctions pour résoudre numériquement l’équation de diffusion :

- `sol_ex(x, y, Lx, Ly)` : calcul de la solution exacte stationnaire.
- Schéma **explicite** pour l’évolution de la vitesse.
- Schéma **implicite** avec matrices creuses (optimisé) et matrices pleines (moins efficace).
- Visualisation 3D de la vitesse à différents instants.

Le script compare les solutions numériques aux solutions exactes pour évaluer la précision des schémas.



## Fonctionnalités

- Calcul de la vitesse $u(x,y,t)$ dans la conduite pour un temps final donné.
- Implémentation de deux schémas numériques :
  - Schéma explicite (conditionnellement stable)
  - Schéma implicite (stable, optimisé avec matrices creuses)
- Application des conditions aux limites (vitesse nulle sur les bords).
- Visualisation graphique 3D des solutions exactes et numériques.
- Calcul et comparaison de l’erreur relative maximale.



## Paramètres du problème

- Dimensions de la conduite : $Lx = 0.02$ m, $Ly = 0.01$ m
- Viscosité cinématique : $\nu = 1 \times 10^{-6}$
- Débit imposé : $Q = 1 \times 10^{-6}$ m³/s
- Discrétisation :
  - $Nx = 16$, $Ny = 8$
  - Pas de temps : $\Delta t = 0.2$
  - Temps final : $tf = 120$ s



## Outils utilisés

- Python
- Numpy
- Scipy
- Matplotlib



## Résultats

- **Solution exacte** : vitesse maximale au centre de la conduite, nulle sur les bords.
- **Schéma explicite** : très proche de la solution exacte, erreur relative maximale $\approx 2.16\%$.
- **Schéma implicite (matrice creuse)** : solution stable et précise, erreur maximale $\approx 1.06\%$.
- **Schéma implicite (matrice pleine)** : mêmes résultats mais temps de calcul très long (~26 minutes), démontrant l’intérêt des matrices creuses.



## Conclusion

Ce TP permet de comparer l’efficacité et la précision de schémas explicites et implicites pour résoudre une équation de diffusion en 2D.  

- Le schéma explicite est simple mais nécessite un petit pas de temps pour rester stable.
- Le schéma implicite est stable pour tout pas de temps, et plus précis lorsqu’il utilise des matrices creuses.
- L’utilisation de matrices creuses optimise fortement le temps de calcul pour des systèmes linéaires de grande taille.



## Auteur

Mouhalhal Moussaab



