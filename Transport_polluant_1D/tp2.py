import numpy as np
import matplotlib.pyplot as plt
#============================================================================================
def thomas(A,D):
    '''
    Calcul de la solution du système linéaire Cx=D dans le cas d'un système tridiagonal en utilisant l'algorithme de thomas (voir wikipedia)
    variables d'entrée:
        A: matrice tridiagonale C mise sous la forme de 3 colonnes et N lignes
        D: vecteur de N composantes
    variable de sortie:
        u: vecteur solution de Au=D
    '''
    N=len(D)
    Dp=np.zeros((N))
    Cp=np.zeros((N))
    u = np.zeros((N))
    Cp[0]=A[0,2]/A[0,1]
    Dp[0]=D[0]/A[0,1]
    for i in range(1,N):
        Cp[i]=A[i,2]/(A[i,1]-A[i,0]*Cp[i-1])
        Dp[i]=(D[i]-A[i,0]*Dp[i-1])/(A[i,1]-A[i,0]*Cp[i-1])
    u[-1]=Dp[-1]
    for i in np.arange(N-2,-1,-1):
        u[i] = Dp[i]-Cp[i]*u[i+1]
    return u
#============================================================================================
def Courant(L,V,beta,t,ntime,N):
    '''
    Calcul de C(x,t) solution de l'équation de transport par le schéma de courant
    '''
    dx= 4*L/ (N-1)           #pas d'espace
    dt=t / ntime      #pas de temps
    r= V*dt/dx  #CFL
    N=int(N)
    ntime=int(ntime)
    #Matrice A mise sous la forme de 3 colonnes et N lignes
    A = np.zeros((N,3))
    A[:,0]=0
    A[:,1]=1
    A[:,2]=0

    #Matrice B mise sous la forme de 3 colonnes et N lignes
    B = np.zeros((N,3))
    B[:,0]=r
    B[:,1]=1- r
    B[:,2]=0

    #Application des conditions aux limites 

    B[0,:]=0
    
    C=np.zeros(N)  #solution initiale
    for i in range(N) :
        x=-L +i*dx
        if abs(x)< L/beta :
            C[i]=np.exp(-1/(1-(beta*x/L)**2))
    I=C #solution initiale
    D = np.zeros(N)
    for k in range(ntime):   #Boucle en temps
        for i in range(1,N):  #On définit le second membre sous la forme d'un vecteur en effectuant B*Un
            D[i]=B[i,0]*C[i-1] + B[i,1]*C[i] 
        D[0]=B[0,1]*C[0] 
        C=thomas(A,D)
    return  C
#============================================================================================
def LaxWendroff(L,V,beta,t,ntime,N):
    '''
    Calcul de C(x,t) solution de l'équation de transport par le schéma de Lax Wendroff
    '''
    dx= 4*L/ (N-1)           #pas d'espace
    dt=t / ntime      #pas de temps
    r= V*dt/dx  #CFL
    N=int(N)
    ntime=int(ntime)
    #Matrice A mise sous la forme de 3 colonnes et N lignes
    A = np.zeros((N,3))
    A[:,0]=0
    A[:,1]=1
    A[:,2]=0

    #Matrice B mise sous la forme de 3 colonnes et N lignes
    B = np.zeros((N,3))
    B[:,0]=0.5*(r+r**2)
    B[:,1]=1-r**2
    B[:,2]=-0.5*(r-r**2)

    #Application des conditions aux limites 
    B[0,:]=0
    B[-1,2]=0
    
    C=np.zeros(N)  #solution initiale
    for i in range(N) :
        x=-L +i*dx
        if abs(x)< L/beta :
            C[i]=np.exp(-1/(1-(beta*x/L)**2))
    I=C #solution initiale
    D = np.zeros(N)
    for k in range(ntime):   #Boucle en temps
        for i in range(1,N-1):  #On définit le second membre sous la forme d'un vecteur en effectuant B*Un
            D[i]=B[i,0]*C[i-1] + B[i,1]*C[i] +B[i,2]*C[i+1]
        D[0]=B[0,1]*C[0] +B[0,2]*C[1]
        D[-1]= B[-1,0]*C[-2] + B[-1,1]*C[-2]
        C=thomas(A,D)
    
    return  C
#============================================================================================
def LeapFrog(L,V,beta,t,ntime,N):
    '''
    Calcul de C(x,t) solution de l'équation de transport par le schéma de Leap Frog
    '''
    dx= 4*L/ (N-1)           #pas d'espace
    dt=t / ntime      #pas de temps
    r= V*dt/dx  #CFL
    N=int(N)
    ntime=int(ntime)
    #Matrice A mise sous la forme de 3 colonnes et N lignes
    A = np.zeros((N,3))
    A[:,0]=0
    A[:,1]=1
    A[:,2]=0

    #Matrice B1 mise sous la forme de 3 colonnes et N lignes
    B1 = np.zeros((N,3))
    B1[:,0]=0
    B1[:,1]=1
    B1[:,2]=0

    #Application des conditions aux limites B1
    B1[0,:]=0

    #Matrice B2 mise sous la forme de 3 colonnes et N lignes
    B2 = np.zeros((N,3))
    B2[:,0]=r
    B1[:,2]=0
    B2[:,2]=-r

    #Application des conditions aux limites B2
    B2[0,:]=0
    B2[-1,2]=0

    
    C0=np.zeros(N)  #solution initiale
    for i in range(N) :
        x=-L +i*dx
        if abs(x)< L/beta :
            C0[i]=np.exp(-1/(1-(beta*x/L)**2))
    I=C0 #solution initiale
    C1=Courant(L,V,beta,dt,ntime,N)
    D = np.zeros(N)
    for k in range(ntime):   #Boucle en temps
        for i in range(1,N-1):  #On définit le second membre sous la forme d'un vecteur en effectuant B*Un
            D[i]=B1[i,1]*C0[i] + B2[i-1,0]*C1[i-1] + B2[i+1,2]*C1[i+1] 
        D[0]=B1[0,1]*C0[0]+ B2[1,2]*C1[1]
        D[-1]=B1[-1,1]*C0[-1] + B2[-2,0]*C1[-2] 
        
        C0=C1
        C1=thomas(A,D)
    return C1
 
#============================================================================================
def sol_ex(x, t, V, L, beta): #Solution analytique

    '''
    Calcul de C(x,t) solution exacte de l'équation de transport
    '''
    C_exact = np.zeros(len(x))
    for i in range(len(x)):
        if abs(x[i] - V * t) < L / beta :
            C_exact[i] = np.exp(-1 / (1 - (beta * (x[i] - V * t) / L)**2))
    return C_exact

#=======================Programme principal===================================================
if __name__ == "__main__":
    #Paramètres physiques:
    L=1
    V=0.3
    beta=5

    #Paramètres liés à la résolution numérique:
    N=500
    dx=4*L/(N-1)
    x=np.linspace(-L ,3*L, N)
        
    #solution numérique:
    
    #Shéma de Courant
    for t in [1, 3, 5] :
        dt=0.5*dx/V     #on impose r=0.5
        ntime=t/dt 
        C=Courant(L,V,beta,t,ntime,N)
        a,=plt.plot(x,C,label='sol num pour t={}'.format(t))
        I=sol_ex(x, t, V, L, beta)
        plt.plot(x,I,'--',color=a.get_color(),markersize=2, label='sol exacte t={} '.format(t))
    plt.title('Shéma de Courant (solution numérique)')
    plt.legend()
    plt.grid(True)
    plt.xlim(-0.3, 2) 
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('courant1.png')
    plt.figure()

    #Shéma de LaxWendroff
    for t in [1, 3, 5] :
        dt=0.5*dx/V
        ntime=t/dt 
        C=LaxWendroff(L,V,beta,t,ntime,N)
        a,=plt.plot(x,C,label='sol num pour t={}'.format(t))
        I=sol_ex(x, t, V, L, beta)
        plt.plot(x,I,'--',color=a.get_color(),markersize=2, label='sol exacte t={} '.format(t))
    plt.title('Shéma de LaxWendroff (solution numérique)')
    plt.legend()
    plt.grid(True)
    plt.xlim(-0.3, 2)  
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('lax1.png')
    plt.figure()


    #Shéma de LeapFrog
    for t in [1, 3, 5] :
        dt=0.5*dx/V
        ntime=t/dt 
        C=LeapFrog(L,V,beta,t,ntime,N)
        a,=plt.plot(x,C,label='sol num pour t={}'.format(t))
        I=sol_ex(x, t, V, L, beta)
        plt.plot(x,I,'--',color=a.get_color(),markersize=2, label='sol exacte t={} '.format(t))
    plt.title('Shéma de LeapFrog (solution numérique)')
    plt.legend()
    plt.grid(True)
    plt.xlim(-0.3, 2)  
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('leapfrog1.png')
    plt.figure()

    #on trace les Schémas pour différentes valeurs de r (t=5s )
    
    #Shéma de Courant
    t=5
    r_=[0.1, 0.5, 0.9]
    for r in r_ :
        dt=r*dx/V
        ntime=t/dt
        C=Courant(L,V,beta,t,ntime,N)
        plt.plot(x,C, label='r={}'.format(round(r,2)))
    I=sol_ex(x, t, V, L, beta)
    plt.plot(x, I,'--', label='sol exacte')
    plt.legend()
    plt.grid(True)
    plt.title('Shéma de Courant')
    plt.xlim(1, 1.8)
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('courant2.png')
    plt.figure()

    #Shéma de LaxWendroff
    for r in r_ :
        dt=r*dx/V
        ntime=t/dt
        C=LaxWendroff(L,V,beta,t,ntime,N)
        plt.plot(x,C, label='r={}'.format(round(r,2)))
    plt.plot(x, I,'--', label='sol exacte')
    plt.legend()
    plt.grid(True)
    plt.title('Shéma de LaxWendroff')
    plt.xlim(1, 1.8) 
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('lax2.png')
    plt.figure()

    #Shéma de LeapFrog
    for r in r_ :
        dt=r*dx/V
        ntime=t/dt
        C=LeapFrog(L,V,beta,t,ntime,N)
        plt.plot(x,C, label='r={}'.format(round(r,2)))
    plt.plot(x, I,'--', label='sol exacte')
    plt.legend()
    plt.grid(True)
    plt.title('Shéma de LeapFrog')
    plt.xlim(1, 1.8) 
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    plt.savefig('leapfrog2.png')
    plt.figure()
        



    
    
