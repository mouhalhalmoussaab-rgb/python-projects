# -*- coding: utf-8 -*-
# auteur NOM:Mouhalhal  PRENOM:Moussaab   NUMERO_ETUDIANT:12408113
# canvas programme Python
# bibliotheque
import numpy as np
import matplotlib.pyplot as plt
# consigne: 
#  - programmation fonction 
#    bien spécifier les arguments et la valeur de retour
#  - verification
#    a la fin du fichier dans la section :
#         if __name__ == '__main__'
#    on appelle la fonction avec un jeux d'arguments pour vérifier

# implementation
def lecture(fichier):
    '''lit les données dans fichier et renvoie les tableaux T,X,Y'''
    f=np.loadtxt(fichier)
    T=f[:,0]
    X=f[:,1]
    Y=f[:,2]
    return T, X, Y

def tri_bulles(T,X,Y):
    '''tri les données sur place (par rapport a T): les tableaux T,X,Y sont donc modifiés'''
    for j in range (len(T)-1 , 1, -1) :
        for i in range(j) :
            if T[i] > T[i+1] :
                T[i],T[i+1] = T[i+1],T[i]
                X[i],X[i+1] = X[i+1],X[i]
                Y[i],Y[i+1] = Y[i+1],Y[i]
    return

def spectre(T,U):
    '''calcul le spectre de U echantillonné en T'''
    dt=T[1]-T[0]
    N=len(U)
    s = np.fft.rfft(U)       
    f = np.fft.rfftfreq(N, dt)
    return f,s

def analyse_spectral(T,U):
    '''calcul mode fondamentale du signal U échantillonné en T '''
    f,s=spectre(T,U)
    i=np.argmax(np.abs(s))
    return np.abs(s[i])/(len(s)-1) , np.angle(s[i]), f[i]

def ecriture(fichier,T,X,Y):
    ''' Ecrire dans fichier de T,X,Y avec le format définit'''
    tri_bulles(T,X,Y)
    f=open(fichier,'w')
    f.write(f"# courbe lissajous N={len(T)} \n")
    for i in range(len(T)) :
        f.write(f'{T[i]} {X[i]} {Y[i]} \n')
    f.close()
    return

    
# test de vérification
if __name__ == '__main__':
    print("vérification lecture")
    T,X,Y=lecture("mon_fichier.dat")
    print(f'T={T[:5]}...')
    print("vérification tri_bulles")
    tri_bulles(T,X,Y)
    check=True
    for i in range(len(T)-1) :
        if T[i]>T[i+1] :
            check=False
    print(check)
    print(f"T trié ={T[:5]} ...")
    print("vérification analyse_spectral")
    t=np.linspace(0,2*np.pi,100, endpoint=False)
    U=2*np.cos(t)
    print(f'Amplitude={analyse_spectral(t,U)[0]}, phase ={analyse_spectral(t,U)[1]}, fréquence={analyse_spectral(t,U)[2]} Hz')
    print("vérification ecriture")
    ecriture("test.dat",T,X,Y)
    f=open("test.dat",'r')
    lines=f.readlines()[:3]
    print(f"{lines} ...")
    

