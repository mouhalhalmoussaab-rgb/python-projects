"""
Projet universitaire - Simulation aérodynamique d’une éolienne à 3 pales (profil NACA4409)
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from xfoil import XFoil
from xfoil.model import Airfoil

#données
nu=15.6*10**-6
rho=1.2
lamda=8
B=3      #nombre de pales
N=20      
V_vent=12
beta=2     #en degré
R=30

#la fonction de la distribution de la corde
a= (0.03-0.1)/(1-0.3)
b= (0.1 - a*0.3)*R
def C(r) :
    C=a*r + b
    return C


#définition des fonction permettant le calcul de la puissance 
omega=V_vent * lamda /R
def U(r) :
    U=r*omega
    return U

V=2/3*V_vent
def phi(r) :
    phi=math.atan2(V, U(r))   #en radian
    phi=phi*180/np.pi         #en degré
    return phi

def alpha(r) :
    alpha=phi(r)-beta
    return alpha

def W(r) :
    W=np.sqrt(U(r)**2 + V**2)
    return W

def Re(r) :
    Re=W(r)*C(r) /nu
    return Re



xf = XFoil()
xf.M = 0.0    # nombre de Mach (0 = subsonique)
xf.iter = 300    # nombre maximal d'itérations pour la convergence
xf.naca(4409)
airfoil = xf.airfoil
x, y = airfoil.x, airfoil.y    #les coordonnées du profil
def cl_cd(xf,Re,alpha):
    xf.Re = Re    # nombre de Reynolds
    xf.print=False
    cl, cd, cm, _ = xf.a(alpha)
    return cl,cd


#la fonction de la puissance élémentaire dP
dr=(R-0.3*R)/(N-1)
def dp(r,cl,cd) :
    phi_=phi(r)*np.pi/180     #en radian
    dp=B*omega*r*(cl*np.sin(phi_) -cd*np.cos(phi_))*1/2 *rho *W(r)**2 *C(r) *dr
    return dp

#récupération de cl et cd 
cl=np.zeros(N)
cd=np.zeros(N)
r=np.linspace(0.3*R, R, N)
for i in range(N) :
    cl[i], cd[i]= cl_cd(xf,Re(r[i]), alpha(r[i]))

#nettoyage de cl et cd par la méthode de la moyenne
for i in range(N) :
    if math.isnan(cl[i])==True :
        cl[i]= (cl[i-1] + cl[i+1])/2 
    if math.isnan(cd[i])==True :
        cd[i]= (cd[i-1] + cd[i+1])/2 



#calcul de la puissance et Cp
P=0
for i in range(len(r)) :
    P+=dp(r[i],cl[i],cd[i])
print('La puissance de l’éolienne pour le profil NACA4409 ={} W'.format( P))
Cp=P /(1/2*rho*V_vent**3 *np.pi *R**2)
print('Cp=',Cp)

#calcul et tracé de alpha, Re et Cp en fonction de r
alpha_=np.zeros(len(r))
Re_=np.zeros(len(r))
cp=np.zeros(len(r))
for j in range(len(r)) :
    Re_[j]=Re(r[j])
    alpha_[j]=alpha(r[j])
    P=0
    for i in range(j) :
        P+=dp(r[i],cl[i],cd[i])
    cp[j]=P /(1/2*rho*V_vent**3 *np.pi *R**2)
m=max(Re_)
i=np.argmax(Re_)
print('le max de Re =', m ) 
print('la valeur de r corresepandente =', r[i] ) 

plt.figure()
plt.plot(r,alpha_,'-o')
plt.title(r'$\alpha$ en fonction du rayon r ')
plt.xlabel('r (m)')
plt.ylabel(r'$\alpha$ (°)')
plt.grid(True)
plt.savefig('alpha.png')

plt.figure()
plt.plot(r,Re_,'-o')
plt.title('Re en fonction du rayon r ')
plt.xlabel('r (m)')
plt.ylabel('Re')
plt.grid(True)
plt.savefig('re.png')

plt.figure()
plt.plot(r,cp,'-o')
plt.title('Cp en fonction du rayon r ')
plt.xlabel('r (m)')
plt.ylabel('Cp')
plt.grid(True)
plt.savefig('cp.png')



#récupération de cl et cd pour 11 valeurs de Re entre 10^96 et 20 . 10^6 pour alpha allant de 0° à 22°
Re=np.linspace(10**6, 20*10**6, 11)
alpha=np.linspace(0,22,23)
NRe=len(Re)
CL=np.zeros((23,NRe))
CD=np.zeros((23,NRe))
for i in range(NRe) :
    xf.Re=Re[i]
    CL[:,i],CD[:,i]=xf.aseq(0, 23,1)[1:3]   #alpha de 0 à 22 degré, 23 valeurs


#nettoyage de la base de données
for j in range(NRe) :
    k=1  #paramètre qui définie si la base est nettoyée (k=0) ou non (k#0)
    while k!=0 :
        #méthode de la moyenne
        for i in range(1,23-1) :   #23 est le nombre de valeurs de alpha 
            if math.isnan(CL[i,j])==True :
                if math.isnan(CL[i-1,j])==False and math.isnan(CL[i+1,j])==False :
                    CL[i,j]= (CL[i-1,j] + CL[i+1,j])/2 

        #approximation linéaire avec les points de droite 
        for i in range(0, 23-2) : 
            if math.isnan(CL[i,j])==True :
                if math.isnan(CL[i+1,j])==False and math.isnan(CL[i+2,j])==False :
                    a=CL[i+2,j]-CL[i+1,j]
                    b=CL[i+1,j]-a*(i+1)
                    CL[i,j]=a*i+b
        #approximation linéaire avec les points de gauche
        for i in range(2, 23) :
            if math.isnan(CL[i,j])==True :
                if math.isnan(CL[i-1,j])==False and math.isnan(CL[i-2,j])==False :
                    a=CL[i-1,j]-CL[i-2,j]
                    b=CL[i-1,j]-a*(i-1)
                    CL[i,j]=a*i+b
        #vérification si si la base est nettoyée ou non
        k=0
        for i in range(0, 23) :
            if math.isnan(CL[i,j])==True :
                k=k+1

for j in range(NRe):
    k = 1
    while k != 0:
        for i in range(1, 23 - 1):
            if math.isnan(CD[i, j]) == True:
                if math.isnan(CD[i - 1, j]) == False and math.isnan(CD[i + 1, j]) == False:
                    CD[i, j] = (CD[i - 1, j] + CD[i + 1, j]) / 2

        for i in range(0, 23 - 2):
            if math.isnan(CD[i, j]) == True:
                if math.isnan(CD[i + 1, j]) == False and math.isnan(CD[i + 2, j]) == False:
                    a = CD[i + 2, j] - CD[i + 1, j]
                    b = CD[i + 1, j] - a * (i + 1)
                    CD[i, j] = a * i + b

        for i in range(2, 23):
            if math.isnan(CD[i, j]) == True:
                if math.isnan(CD[i - 1, j]) == False and math.isnan(CD[i - 2, j]) == False:
                    a = CD[i - 1, j] - CD[i - 2, j]
                    b = CD[i - 1, j] - a * (i - 1)
                    CD[i, j] = a * i + b
        k = 0
        for i in range(0, 23):
            if math.isnan(CD[i, j]) == True:
                k = k + 1

#tracé du profil
plt.figure()
plt.plot(x,y)
plt.title('Profil naca4409')
plt.axis('equal')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('profile.png')

n=1
#tracé de CL en fonction de alpha pour Re fixe = Re[n]
plt.figure()
plt.plot(alpha,CL[:,n],'-o')
plt.title(r'CL en fonction $\alpha$ pour Re={}'.format(Re[n]))
plt.xlabel(r'$\alpha$ (°)')
plt.ylabel('CL')
plt.grid(True)
plt.savefig('cl.png')
#max de cl 
m=max(CL[:,n])
i=np.argmax(CL[:, n])
print('le max de Cl =', m ) 
print('la valeur de alpha corresepandente =', alpha[i] ) 


#tracé de CD en fonction de alpha pour Re fixe = Re[n]
plt.figure()
plt.plot(alpha,CD[:,n],'-o')
plt.title(r'CD en fonction $\alpha$ pour Re={}'.format(Re[n]))
plt.xlabel(r'$\alpha$ (°)')
plt.ylabel('CD')
plt.grid(True)
plt.savefig('cd.png')
