from math import *
from functions import *
'''

def hist_pos(N):
    m=5
    simulation_RPMD(N,Vp,Vpprime)
    X=np.arange(-m*a,m*a,e)
    print(len(X))
    Y=np.zeros((len(X)))
    for i in range(pas):

        for j in range(N):
            x=0
            x+=atome[j].x[i]
            Y[int((x+m*a)/e)]+=1/N
    return X,Y

'''
def hist_pos(N):
    m=5
    simulation_RPMD(N,Vp,Vpprime)
    X=np.arange(-m*a,m*a,e)
    Y=np.zeros((len(X)))
    for i in range(pas):

        for j in range(N):
            x=0
            x+=atome[j].x[i]
            Y[int((x+m*a)/e)]+=1/N
    return X,Y




X,Y=hist_pos(P)
plt.plot(X,Y)
plt.xlabel("abscisse (en m)")
plt.ylabel("Nombre de points")
plt.title("Histogramme des positions prises par les particules pour P="+str(P))
plt.show()
