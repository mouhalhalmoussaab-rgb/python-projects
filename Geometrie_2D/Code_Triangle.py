#pour tester : test_exo Triangle.py exo510 Triangle Point
import numpy as np
from Point import *
# classe Triangle
class Triangle(object):
    """creation d'un triangle dans l'espace cartesien"""
    def __init__(self,P1,P2,P3):
        ''' initialisation '''
        # attention on fait une copie des points !!!
        self.Pts = [P1.copy(),P2.copy(),P3.copy()]
        return
    def __str__(self):
        '''conversion chaine pour affichage'''
        return "Triangle:(%s,%s,%s)"%(self.Pts[0],self.Pts[1],self.Pts[2]) 
    def __eq__(self,T):
        '''test si le triangle est confondu avec T'''
        check=0
        for i in range(3) :
            for j in range(3) :
                if self.Pts[i]==T.Pts[j] :
                    check+=1
        if check==3 :
            return True
        return False
    def barycentre(self):
        '''calcul le point barycentre du triangle'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        x=(p1.x+p2.x+p3.x)/3
        y=(p1.y+p2.y+p3.y)/3
        return Point(x,y)
    def perimetre(self):
        '''calcul perimetre du triangle'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        return p1.distance(p2)+p2.distance(p3)+p3.distance(p1)
    def surface(self):
        '''calcul la surface du triangle'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        return abs(produit_vect(p1,p2,p3))/2
    def coordbary(self,P):
        '''les coordonnées barycentriques'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        
        b1=produit_vect(P,p2,p3)/2/self.surface()
        b2=produit_vect(P,p3,p1)/2/self.surface()
        b3=1-b1-b2
        return b1,b2,b3
    def inside(self,P):
        '''test si le point P est à l intérieur du triangle'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        s=Triangle(P,p2,p3).surface() + Triangle(P,p3,p1).surface() +Triangle(P,p1,p2).surface()
        if s==self.surface() :
            return True
        return False
    def translation(self,dx,dy):
        '''translation du triangle de dx,dy'''
        self.Pts[0].translation(dx,dy)
        self.Pts[1].translation(dx,dy)
        self.Pts[2].translation(dx,dy)
        return 
    def rotation(self,O,alpha):
        """rotation du triangle de alpha autour de O"""
        self.Pts[0].rotation(O,alpha)
        self.Pts[1].rotation(O,alpha)
        self.Pts[2].rotation(O,alpha)
        return
    def rayon(self):
        '''calcul le rayon du cercle contenant le triangle et centré au barycentre'''
        p1=self.Pts[0]
        p2=self.Pts[1]
        p3=self.Pts[2]
        b=self.barycentre()
        return np.max([distance(b,p1), distance(b,p2), distance(b,p3)])
        
    def intersection(self,T):
        '''détermine si il y a une intersection non nulle avec le traingle T'''
        p=self.Pts
        pt=T.Pts
        indices=[0,1,2,0]
        for i in range(3) :
            n1,n2=indices[i],indices[i+1]
            for j in range(3) :
                n3,n4=indices[j],indices[j+1]
                #domaine d'intersection
                x_min,x_max=xmin_xmax(p[n1],p[n2],pt[n3],pt[n4])
                y_min,y_max=ymin_ymax(p[n1],p[n2],pt[n3],pt[n4])
                if x_min<=x_max and y_min<=y_max : #donc il exisite un domaine d'intersection et on vérifiera s'il y a des points d'intersection
                    if p[n1].x==p[n2].x : #si un coté est verticale donc on aura pas l'équation de droite y=a*x+b
                        x=p[n1].x
                        if pt[n3].x==pt[n4].x :
                            if x==pt[n3].x :
                                return True
                        else :
                            a2,b2=droite(pt[n3],pt[n4])
                            y=a2*x+b2
                            if y>=y_min and y<=y_max :
                                return True
                    else :
                        a1,b1=droite(p[n1],p[n2])
                        if pt[n3].x==pt[n4].x :
                            x=pt[n3].x
                            y=a1*x+b1
                            if y>=y_min and y<=y_max :
                                return True
                        else :
                            a2,b2=droite(pt[n3],pt[n4])
                            if a2==a1 :
                                if b2==b1 :
                                    return True
                            else :
                                x=- (b2-b1)/(a2-a1)
                                if x>=x_min and x<=x_max :
                                    return True
        return False
    def plot(self,col='g'):
        '''tracer du triangle avec remplissage'''
        p=self.Pts
        x = [p[0].x, p[1].x, p[2].x, p[0].x]
        y = [p[0].y, p[1].y, p[2].y, p[0].y]
        plt.fill(x, y, color=col, alpha=0.3, label=f"{self}")
        plt.legend()
        plt.grid(True)
        plt.title("Tracé des triangles")
        return

def droite(p1,p2) :
    '''les coefficients a et b tel que y=a*x+b'''
    a= (p2.y-p1.y) / (p2.x-p1.x)
    b=p1.y-a*p1.x
    return a,b
    
def xmin_xmax(p1,p2,p3,p4) :
    #intervalle de x du coté p1p2 de T1
    x_min1=min(p1.x,p2.x)
    x_max1=max(p1.x,p2.x)
    #intervalle de x du coté p3p4 de T2
    x_min2=min(p3.x,p4.x)
    x_max2=max(p3.x,p4.x)
    #le max des min 
    x_min=max(x_min1,x_min2)
    #le min des max
    x_max=min(x_max1,x_max2)
    return x_min, x_max

def ymin_ymax(p1,p2,p3,p4) :
    y_min1=min(p1.y,p2.y)
    y_max1=max(p1.y,p2.y)
    y_min2=min(p3.y,p4.y)
    y_max2=max(p3.y,p4.y)
    y_min=max(y_min1,y_min2)
    y_max=min(y_max1,y_max2)
    return y_min, y_max
    
# test de vérification
if __name__ == '__main__':
    A=Point(0,1)
    B=Point(1,0)
    C=Point(0,0)
    D=Point(1,2)
    E=Point(1.9,1)
    F=Point(1,1)
    T=Triangle(Point(0.0,0.0),Point(1.0,0.0),Point(0.0,1.0))
    T1=Triangle(Point(1,0.5),Point(-1,0.5),Point(-1,1.0))
    print("vérification de la méthode intersection : ",T.intersection(T1))
    T.plot()
    T1.plot(col='r')
