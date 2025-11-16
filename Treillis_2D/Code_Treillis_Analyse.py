from treillis import *

# Lecture des données du treillis à partir du fichier "treillis.dat"
tr = lecture("treillis.dat")

# Propriétés matérielles du treillis
tr.E = 200*1e9       # Module de Young en Pa
tr.rho = 8000        # Masse volumique en kg/m³

# Géométrie des barres
d = 8*1e-2           # Diamètre extérieur des barres en m
e = 1e-2             # Épaisseur des barres en m (barres trouées)
tr.S = np.pi*(d**2-e**2)/4  # Section transversale des barres

# Paramètres de charge
ec = 4               # Écartement entre deux treillis similaires en m
k = 80               # Poids de la toiture en kg/m²
hauteur_neige = 11e-2  # Hauteur de neige en m
rho_neige = 300        # Masse volumique de la neige en kg/m³
g = 9.81              # Accélération gravitationnelle en m/s²

# Affichage du treillis initial
trace_treillis(tr,'Treillis initial')

# Calcul de la matrice de rigidité globale
AA = matrice_rigidite(tr)

# Définition des conditions aux limites
tr.CL = np.zeros(2*tr.nn)  # Initialisation des CL (0 = libre)
tr.CL[0] = 3               # Nœud 0 encastré (ux=uy=0)
tr.CL[8] = 3               # Nœud 8 encastré
tr.FCL = np.zeros((tr.nn,2))  # Initialisation des forces nodales

# Application du poids propre des barres
for i in range(tr.ne):
    n1, n2 = tr.G[i]
    poids = -tr.S * long_barres(tr)[i] * tr.rho * g  # Force verticale négative
    tr.FCL[n1][1] += poids/2  # Répartition moitié/moitié sur les nœuds
    tr.FCL[n2][1] += poids/2

# Application du poids de la toiture et de la neige sur certaines barres
for i in [5, 4, 3, 13, 14, 15]:
    n1, n2 = tr.G[i]
    poids_toiture = -long_barres(tr)[i] * ec * k * g      # Poids toiture par barre
    poids_neige = -hauteur_neige * long_barres(tr)[i] * ec * rho_neige * g  # Poids neige par barre
    # Répartition sur les quatre nœuds de la surface créée par la barre et son symétrique
    tr.FCL[n1][1] += poids_toiture/4 + poids_neige/4
    tr.FCL[n2][1] += poids_toiture/4 + poids_neige/4

# Application des conditions aux limites sur la matrice et le vecteur des forces
A, B = conditions_limites(tr, AA)

# Résolution du système linéaire pour obtenir les déplacements
U = np.linalg.solve(A, B)

# Affichage du déplacement maximum
print(f"Le max du déplacement= { np.max(np.abs(U))} m")

# Tracé du treillis sous effet des charges (déplacements amplifiés x500)
trace_treillis(tr,"Treillis sous l'effet des charges appliquées (déplacements x500)", U*500)
plt.savefig("treillis.pdf")  

# Calcul des nouvelles coordonnées des nœuds après déplacement
C = tr.X + U.reshape(np.shape(tr.X))

# Détermination de l'état de chaque barre (traction ou compression)
for i in range(len(tr.G)):
    n1, n2 = tr.G[i]
    x = C[n2][0] - C[n1][0]
    y = C[n2][1] - C[n1][1]
    norme = np.sqrt(x**2 + y**2)
    if norme < long_barres(tr)[i]:
        print("la barre {} est en compression".format(i))
    else:
        print("la barre {} est en Traction".format(i))
