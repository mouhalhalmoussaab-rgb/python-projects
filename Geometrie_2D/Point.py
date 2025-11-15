import numpy as np
import matplotlib.pyplot as plt
# precision des calcul geometrique
EPS_GEO = 1.e-5
# classe Point
class Point(object):
    """création d'un point dans l'espace cartesien"""
    def __init__(self,x,y):
        ''' initialisation '''
        self.x = x
        self.y = y
        return 
    def copy(self):
        '''renvoie une copie du point'''
        return Point(self.x,self.y)
    def __str__(self):
        '''conversion chaine pour affichage'''
        return "Point:(%s,%s)"%(self.x,self.y) 
    def __eq__(self,P):
        ''' test si le point est confondu avec P a epsilon pres'''
        return distance(self,P) < EPS_GEO
    def rayon(self):
        '''calcul le rayon du point / origine'''
        return distance(self,Point(0,0))
    def angle(self):
        '''calcul angle en degré de -180 a 180'''
        return np.arctan2(self.y,self.x) * 180/np.pi
    def distance(self,P):
        '''calcul la distance au point P'''
        return distance(self,P)
    def rotation(self, O, beta):
        '''Rotation du point par rapport au point O avec angle beta (en degré)'''
        beta *=np.pi/180   # conversion en radians
        B = np.array([self.x - O.x, self.y - O.y])
        # matrice de rotation
        R = np.array([[np.cos(beta), -np.sin(beta)],
                        [np.sin(beta),  np.cos(beta)]])
        A = R @ B
        self.x = O.x + A[0]
        self.y = O.y + A[1]
    def translation(self,dx,dy):
        '''translation du point de dx,dy'''
        self.x+=dx
        self.y+=dy
        return
    def polaire(self):
        '''coordonnées polaire r,theta (en radian) du point'''
        return self.rayon(), self.angle()*np.pi/180
    def plot(self):
        '''tracer du point'''
        plt.plot(self.x,self.y,'o',label=f"{self}")
        plt.legend()
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("tracé des Points")
        plt.grid(True)
        return
#
# Fonctions utiles
#
# distance entre 2 points
def distance(P1,P2):
    """ calcul distance entre 2 points"""
    return np.sqrt( (P1.x-P2.x)**2 + (P1.y-P2.y)**2 )
# produit vectoriel AB x AC
def produit_vect(A,B,C):
    """ calcul produit vectoriel AB.AC (cpste suivant z)"""
    return (B.x-A.x)*(C.y-A.y)-(B.y-A.y)*(C.x-A.x)

# test de vérification
if __name__ == '__main__':
    print("vérification __eq__")
    A=Point(0,1)
    B=Point(0,1)
    C=Point(1,0)
    print(A==B)
    print("vérification rayon")
    print(A.rayon())
    print("vérification angle")
    print(A.angle())
    print("vérification distance")
    print(A.distance(B))
    print("vérification rotation")
    A.rotation(C,90)
    print(A)
    print("vérification translation")
    B.translation(1,1)
    print(B)
    print("vérification polaire")
    D=Point(-1,0)
    print(f"rayon= {D.polaire()[0]}, angle ={D.polaire()[1]} rad")
    print("vérification plot")
    E=Point(3,4)
    C.plot()
    D.plot()
    C.plot()
    
