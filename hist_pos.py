from math import *
from functions import *


def hist_pos(N):
    simulation_RPMD(N,Vp,Vpprime)
    X=np.arange(-a,a,e*10)
    print(len(X))
    Y=np.zeros((len(X)))
    for i in range(pas):

        for j in range(N):
            x=0
            x+=atome[j].x[i]
            Y[int((x+m*a)/10/e)]+=1/N
    return X,Y






X,Y=hist_pos(P)
print(Y)
plt.hist(Y,X)
plt.xlabel("abscisse (en m)")
plt.ylabel("Nombre de points")
plt.title("Histogramme des positions prises par les particules pour P="+str(P))
plt.show()
