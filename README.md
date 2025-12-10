# Python Projects

Ce dépôt contient plusieurs projets Python réalisés dans le cadre de mon Master en Simulation et Modélisation Mécanique.  
Chaque projet contient son propre README détaillant la méthodologie, les outils utilisés et les résultats.

## Projets inclus

- [Transport_polluant_1D](./Transport_polluant_1D) – Simulation numérique de l’équation de transport 1D pour un polluant se déplaçant à vitesse constante.  
  Le projet inclut la comparaison de la solution exacte avec trois schémas numériques : Courant, Lax-Wendroff et Leap-Frog. Les scripts calculent la concentration C(x,t), tracent les profils pour différents instants et étudient l’effet du paramètre CFL (r) sur la diffusion et la dispersion numérique.

- [Ecoulement_instationnaire_2D](./Ecoulement_instationnaire_2D) – Étude de la diffusion instationnaire en 2D dans une conduite rectangulaire.  
  Le projet implémente des schémas explicites et implicites pour résoudre l’équation de diffusion. Il compare les solutions numériques aux solutions exactes, visualise les champs de vitesse et illustre l’importance des matrices creuses pour optimiser le temps de calcul dans le cas de grandes matrices.
 

- [Meteo_Analysis_Pandas](./Meteo_Analysis_Pandas) – Analyse des relevés météo d’Ajaccio pour l’année 2014 avec Python et Pandas.  
  Le projet inclut le traitement du DataFrame, l’étude des températures, précipitations, neige et vent, ainsi que la comparaison avec les normales climatiques.
  
- [Geometrie_2D](./Geometrie_2D) – Projet de calcul géométrique en 2D portant sur la manipulation de points et de triangles.  
  Il inclut l’implémentation des opérations géométriques (distances, barycentres, surfaces), des transformations (translation, rotation), des tests d’intersection entre triangles, ainsi que des visualisations et notebooks de validation.

- [Aerodynamique_XFoil](./Aerodynamique_XFoil) – Simulation aérodynamique d’une éolienne à 3 pales (profil NACA4409).  
  Le projet comprend le calcul des coefficients aérodynamiques (CL, CD), la construction d'une base de données propre et l’évaluation des performances (CP, puissance).

- [Treillis_2D](./Treillis_2D) – Étude et analyse de treillis 2D à l’aide de scripts Python. Le projet inclut la modélisation des structures à barres, le calcul des matrices de rigidité, l’application des conditions aux limites et l’analyse des déplacements sous différentes charges (poids propre, toiture, neige). Un treillis simple est utilisé pour valider le code, puis un treillis complet est analysé pour déterminer les sollicitations de chaque barre (traction ou compression).
  
- [Traitement_signaux](./Traitement_signaux) – Analyse et filtrage de signaux périodiques bruités X(t) et Y(t) pour extraire les composantes principales et les coefficients de Lissajous.  
  Le projet inclut la lecture et le tri des données, le calcul des spectres via FFT, la reconstruction des signaux filtrés et la visualisation des courbes de Lissajous avant et après filtrage. 

## Auteur
Mouhalhal Moussaab

