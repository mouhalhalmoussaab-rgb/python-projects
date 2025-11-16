# Treillis 2D

Ce projet contient des scripts Python permettant d'étudier le comportement mécanique d'un treillis 2D sous différentes charges.

## Description

Le projet est organisé en deux scripts principaux :

- **Code_Treillis.py** : définition de la classe `Treillis` et des fonctions principales pour :
  - Lire les données d'un treillis à partir d'un fichier,
  - Construire la matrice de rigidité,
  - Appliquer les conditions aux limites (nœuds libres, bloqués ou encastrés),
  - Calculer les déplacements des nœuds.
  
  Un exemple de treillis simple est inclus pour valider le fonctionnement du code.

- **Code_Treillis_Analyse.py** : étude d'un treillis complet soumis à plusieurs charges appliquées selon l'axe vertical \(y\) :
  - Poids propre des barres,
  - Poids de la toiture,
  - Poids de la neige.
  
  Le script calcule les déplacements et permet de déterminer si chaque barre est en traction ou en compression.

## Fonctionnalités

- Calcul des déplacements des nœuds d’un treillis 2D.
- Visualisation du treillis initial et déformé.
- Détermination du type de sollicitation (traction/compression) pour chaque barre.
- Support pour plusieurs types de conditions aux limites par nœud :
  - Libre, bloqué en \(x\), bloqué en \(y\), encastré.

## Utilisation

1. Placer les fichiers de données du treillis dans le même répertoire que les scripts.
2. Exécuter `Code_Treillis.py` pour valider le code avec un treillis simple.
3. Exécuter `Code_Treillis_Analyse.py` pour analyser un treillis complet et visualiser les résultats.

## Outils utilisés

- Python (Numpy, Matplotlib)

## Résultats

- Déplacement des nœuds du treillis.
- Identification des barres en traction ou en compression.
- Graphiques du treillis initial et déformé.

## Auteur

Mouhalhal Moussaab
