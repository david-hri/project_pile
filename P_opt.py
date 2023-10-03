from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt
X=[]
V=Vp
P=10
Vprime=Vpprime

def Moy(X):
    x=0
    for i in X:
        x+=i
    return x/len(X)

def P_opt(n):
    X=[]
    for i in range(2,n):
        print(i)

        simulation_RPMD(i,V,Vprime)
        A=atome[0].x
        for j in range(1,i):
        
            A+=atome[j].x
        for k in range(len(A)):
            A[k]/=i
        X.append(Veff(A,V))
    return X

def nombre_plus_proche(L,cible=0.95):
    dernier_element = L[-1]
    objectif = cible * dernier_element
    i=0
  
    ecart_minimal = abs(L[0] - objectif)

    for j in range(len(L)):
        ecart = abs(L[j] - objectif)
        if ecart < ecart_minimal:
            ecart_minimal = ecart
            i=j

    return i

X=P_opt(12)
print(X)
print(nombre_plus_proche(X,0.98)+1)
plt.plot([i+1 for i in range(len(X))],X,label="T="+str(T))
plt.legend()
plt.show()