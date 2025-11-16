#! /usr/bin/env python3
# Bibliothèque pour l'étude d'un treillis 2D
# (C) MOUHALHAL MOUSSAAB 12408113  
import numpy as np
import matplotlib.pyplot as plt

# Définition de la structure Treillis
class Treillis():
    """ Structure de données pour un treillis plan 2D """
    def __init__(self, mon_fichier):
        self.fichier = mon_fichier
        # Nombre de nœuds
        self.nn = 0
        # Coordonnées des nœuds (matrice nn x 2)
        self.X  = None
        # Conditions aux limites (code 0,1,2,3)
        self.CL = None
        # Forces nodales suivant x et y (matrice nn x 2)
        self.FCL = None
        # Nombre de barres
        self.ne = 0
        # Numéros des deux nœuds de chaque barre
        self.G  = None
        # Section, module de Young et masse volumique des barres (constantes)
        self.S  = None
        self.E  = None
        self.rho = None

# Fonction pour lire les données d'un treillis depuis un fichier
def lecture(nom_fichier):
    """ Lecture des données d'un treillis depuis le fichier """
    tr = Treillis(nom_fichier)
    f = open(nom_fichier, 'r')
    lines = f.readlines()
    h = lines[1].split()
    ne = int(h[0])  # nombre de barres
    nn = int(h[1])  # nombre de nœuds
    X = np.zeros((nn,2))      # coordonnées des nœuds
    G = np.zeros((ne,2),dtype=int)  # connectivité des barres
    i = 0
    # Lecture des coordonnées des nœuds
    for line in lines[3:3+nn]:
        X[i] = line.split()
        i += 1
    j = 0
    # Lecture des barres (connectivité)
    for line in lines[4+nn:4+nn+ne]:
        G[j] = line.split()
        j += 1
    # Affectation des valeurs au treillis
    tr.nn = nn
    tr.ne = ne
    tr.X = X
    tr.G = G
    return tr

# Calcul de la longueur de chaque barre
def long_barres(tr):
    """ Calcule un vecteur contenant la longueur des barres du treillis """
    h = np.zeros(tr.ne)
    for i in range(len(tr.G)):
        n1, n2 = tr.G[i]
        h[i] = np.sqrt((tr.X[n1][0] - tr.X[n2][0])**2 + (tr.X[n1][1] - tr.X[n2][1])**2)
    return h

# Matrice de rotation pour passer du repère local au repère global
def rotation_barre(tr,l):
    """ Calcule la matrice de rotation 4x4 pour la barre l """
    L = long_barres(tr)[l]
    n1, n2 = tr.G[l]
    adj = tr.X[n2][0] - tr.X[n1][0]  # côté adjacent
    op = tr.X[n2][1] - tr.X[n1][1]   # côté opposé
    cos = adj / L
    sin = op / L
    T = np.array([[ cos,sin, 0,0],
                  [-sin,cos, 0,0],
                  [ 0,0, cos,sin],
                  [ 0,0,-sin,cos]])
    return T

# Assemblage de la matrice de rigidité globale
def matrice_rigidite(tr):
    """ Calcule la matrice de rigidité globale avant application des CL """
    E = tr.E
    S = tr.S
    nn = tr.nn
    G = tr.G
    A = np.zeros((2*nn,2*nn))
    for i in range(len(tr.G)):
        L = long_barres(tr)[i]
        
        # Matrice de rigidité locale
        K = E*S/L * np.array([[1,0,-1,0],
                              [0,0,0,0],
                              [-1,0,1,0],
                              [0,0,0,0]]) 
        
        # Transformation au repère global
        Ae = rotation_barre(tr,i).T @ K @ rotation_barre(tr,i)
        
        # Assemblage dans la matrice globale
        n1, n2 = G[i]
        ddl = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
        for m in range(4):
            for n in range(4):
                I = ddl[m]
                J = ddl[n]
                A[I][J] += Ae[m][n]
    return A

# Application des conditions aux limites
def conditions_limites(tr,AA):
    """ Applique les CL sur la matrice AA et retourne A et B modifiées """
    # cl=0 aucune CL, cl=1 ux=0, cl=2 uy=0, cl=3 encastrement complet
    cl = tr.CL
    fcl = tr.FCL
    nn = tr.nn
    A = AA.copy()
    B = np.zeros(2*nn)
    for i in range(nn):
        if cl[i]==1:  # ux=0
            A[2*i][2*i]=1
            for j in range(2*nn):
                if j != 2*i:
                    A[2*i][j] = 0
                    A[j][2*i] = 0
            B[2*i]=0
            B[2*i+1]=fcl[i][1]
        if cl[i]==2:  # uy=0
            A[2*i+1][2*i+1]=1
            for j in range(2*nn):
                if j != 2*i+1:
                    A[2*i+1][j]=0
                    A[j][2*i+1]=0
            B[2*i]=fcl[i][0]
            B[2*i+1]=0
        if cl[i]==3:  # encastrement complet
            A[2*i][2*i]=1
            A[2*i+1][2*i+1]=1
            for j in range(2*nn):
                if j!=2*i:
                    A[2*i][j]=0
                    A[j][2*i]=0
                if j!=2*i+1:
                    A[2*i+1][j]=0
                    A[j][2*i+1]=0
            B[2*i]=0
            B[2*i+1]=0
        else:  # nœud libre
            B[2*i] = fcl[i][0]
            B[2*i+1] = fcl[i][1]
    return A, B

# Fonction pour tracer le treillis
def trace_treillis(tr, titre, U=None, style='-g'):
    """ Tracer le treillis avec ou sans déplacements """
    Xd = tr.X.copy()
    G = tr.G
    if U is not None:
        style='-r'  # couleur rouge si déformé
        U = U.reshape(np.shape(tr.X))
        Xd[:,0] = Xd[:,0] + U[:,0]
        Xd[:,1] = Xd[:,1] + U[:,1]
        
    if U is None:
        plt.figure()
        # Tracer les nœuds
        plt.plot(Xd[:,0], Xd[:,1],'o',markersize=7,color='b')
        for i in range(tr.ne):
            plt.text((Xd[G[i][0],0]+Xd[G[i][1],0])/2,
                     (Xd[G[i][0],1]+Xd[G[i][1],1])/2,
                     str(i),fontsize=12,color='b')
        for i in range(tr.nn):
            plt.text(Xd[i,0],Xd[i,1],str(i),fontsize=12)

    # Tracer les barres
    for i in range(tr.ne):
        n1 = tr.G[i,0]
        n2 = tr.G[i,1]
        plt.plot([Xd[n1,0],Xd[n2,0]],[Xd[n1,1],Xd[n2,1]],style,lw=3)
    plt.axis('equal')
    plt.title(titre)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.grid(True)
    return

# Validation du code si le script est exécuté directement
if __name__ == "__main__":
    print("Test lecture")
    tr = lecture('treillis_simple.dat')
    tr.S = 0.000025
    tr.E = 200*1e9
    tr.rho = 8000
    print("Treillis {}  \n nbr de nds={}  \n nbr de barres={}  \n coordonnées des noeuds \n {}  \n éléments \n {} ".format(tr.fichier,tr.nn,tr.ne,tr.X,tr.G))

    print("\n\nTest long_barres")
    print(long_barres(tr))

    print("\n\nTest rotation_barre (barre numéro 0)")
    print(rotation_barre(tr,0))
    
    print("\n\nTest matrice_rigidite")
    print(np.round(matrice_rigidite(tr),2))
    
    print("\n\nTest conditions_limites")
    AA = matrice_rigidite(tr)
    tr.CL = np.array([3, 0, 3])
    F = -1000    # force appliquée
    tr.FCL = np.array([[0, 0], [0, F], [0, 0]])  # exemple de forces
    A,B = conditions_limites(tr,AA)  # application des CL
    print("Matrice de rigidité après application des CL \n {} \nSeconde membre après application des CL \n {}".format(np.array2string(A,suppress_small=True,precision=2),B))

    # Résolution du système linéaire
    U = np.linalg.solve(A,B)
    U_rdm = F*long_barres(tr)[0]/(tr.E*tr.S)
    print(f"\n\nComparaison : le déplacement calculé par le script {U[3]} m, la solution théorique de la RDM {U_rdm} m")

    print("\n\nTest trace_treillis")
    trace_treillis(tr,"Tracé du treillis")
    trace_treillis(tr,"Treillis simple sous l'effet du poids (déplacements x500)", U*500)
    plt.savefig("treillis_simple.pdf")
