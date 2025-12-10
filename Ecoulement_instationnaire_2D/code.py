import numpy as np
import scipy
import matplotlib.pyplot as plt
import time
Nx=16
Ny=8
tf=120
Lx=0.02 #m
Ly=0.01 #m
nu=1e-6
Dt=0.2
Nt=tf/Dt
Nt=int(Nt)
Q=1e-6

def sol_ex(x,y,Lx,Ly) :
    s1=0
    s2=0
    for i in range(25) :
        n=2*i+1
        s1+=n**-4 *(1- 1/np.cosh(n/2 *np.pi*Ly/Lx))
        s2+= n**-3 * (1- np.cosh(n*np.pi*(y-Ly/2)/Lx) / np.cosh(n/2*np.pi*Ly/Lx) ) *np.sin(n*np.pi*x/Lx)
    G=np.pi**4 *Q /( 4*Lx**3 *Ly *s1)
    u_ex=4*G*Lx**2 / np.pi**3 * s2
    return u_ex ,G


#trcé de la sol exacte
x=np.linspace(0,Lx, Nx)
y=np.linspace(0,Ly, Ny)
X,Y=np.meshgrid(x,y)

Z=sol_ex(X,Y,Lx,Ly)[0]
#cont=plt.contour(X,Y,Z)
#plt.colorbar(cont)
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
surf=ax.plot_surface(X,Y,Z, cmap='viridis')
fig.colorbar(surf)

Dx=Lx/(Nx-1)
Dy=Ly/(Ny-1)

G=sol_ex(X,Y,Lx,Ly)[1]
cte_x=Dt*nu /Dx**2
cte_y=Dt*nu /Dy**2 
a=cte_y * np.ones(Nx* Ny -Nx)
b=cte_x * np.ones(Nx* Ny -1)
c=(1 - 2*cte_x - 2*cte_y) * np.ones(Nx* Ny )
d=b
e=a

#sol num schéma explicite
B= np.diag(a,-Nx) + np.diag(b,-1) + np.diag(c) + np.diag(d,1) + np.diag(e,Nx)
C= nu*G*Dt * np.ones(Nx*Ny)



U=np.zeros(Nx*Ny)
for k in range(Nt) :      #choix de Nt compare les maxe des 2 si inf à 1 pour cent c'est parfait
    U= B@U + C 
    for j in range(Ny) :
        U[j*Nx]=0
        U[Nx-1 + j*Nx]=0
    for i in range(Nx) :
        U[i]=0
        U[i + (Ny-1)*Nx]=0
U=U.reshape(X.shape)
    

Z=U
#print('err schéma explicite = ', abs(max(Z)-max(sol_ex(X,Y,Lx,Ly)[0]))/max(sol_ex(X,Y,Lx,Ly)[0]))
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
surf=ax.plot_surface(X,Y,Z, cmap='viridis')
fig.colorbar(surf)

#sol num schéma implicite
#t1=time.process_time()
a=-cte_y * np.ones(Nx* Ny -Nx)
b=-cte_x * np.ones(Nx* Ny -1)
c=(1 + 2*cte_x + 2*cte_y) * np.ones(Nx* Ny )
d=b
e=a
A=np.diag(a,-Nx) + np.diag(b,-1) + np.diag(c) + np.diag(d,1) + np.diag(e,Nx)
A=scipy.sparse.csr_matrix(A)

U=np.zeros(Nx*Ny)
for k in range(Nt) :      
    U=scipy.sparse.linalg.spsolve(A,U+C)
    for j in range(Ny) :
        U[j*Nx]=0
        U[Nx-1 + j*Nx]=0
    for i in range(Nx) :
        U[i]=0
        U[i + (Ny-1)*Nx]=0
U=U.reshape(X.shape)

#t2=time.process_time()
#duree=t2-t1
    

Z=U
#print('err schéma implicite = ', abs(max(Z)-max(sol_ex(X,Y,Lx,Ly)[0]))/max(sol_ex(X,Y,Lx,Ly)[0]))
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
surf=ax.plot_surface(X,Y,Z, cmap='viridis')
fig.colorbar(surf)
