Traitement et filtrage de signaux


## Objectif
Analyser et filtrer deux signaux périodiques X(t) et Y(t) représentant des courbes de Lissajous bruitées afin d’extraire les caractéristiques principales (fréquences, amplitudes, coefficients de Lissajous).

## Tâches principales
- Lecture des données depuis un fichier texte (`donnes.dat`)
- Tri des données par rapport au temps
- Calcul des spectres des signaux avec la transformée de Fourier (FFT)
- Extraction des amplitudes, phases et fréquences fondamentales
- Reconstruction des signaux filtrés X1 et Y1
- Visualisation des signaux et des courbes de Lissajous avant et après filtrage

## Outils utilisés
- Python (NumPy, FFT, Matplotlib)
- Méthodes numériques : FFT, tri, calcul de spectres, visualisation graphique

## Résultats
- Spectre de X : signal bruité avec f ≈ 0.125 Hz  
- Spectre de Y : signal peu bruité avec f ≈ 0.143 Hz  
- Signaux filtrés X1 et Y1 plus réguliers
- Coefficients de Lissajous : p ≈ 2.759, q ≈ 1.747

## Auteur
Mouhalhal Moussaab

## Date
21/11/2025
