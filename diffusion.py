from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt
import scipy


n=20        # Nombre de valeurs de t (t=j*dt) pour lesquelles ont calcule la valeur moyenne 
dj=10000    # Ecart entre deux valeurs de j successives pour lesquelle on calcule l'écart quadratique moyen
c=1         # Permet de décaler la première valeur de j (j_0=dj*c) pour laquelle on calcule l'écart quadratique moyen


simulation_RPMD(P,Vsinprime)

def moyenne_decalage_de_l_origine (j) :   # Calcule l'écart quadratique moyen pour j fixé : moy((x(t_0+t)-x(t_0)**2) avec t_0=k*dt et t=j*dt
    S = 0
    for k in range (pas-j):
        S+=(atome[0].x[j+k]-atome[0].x[k])**2
    Moy_dec = S/(pas-j)
    return Moy_dec

def  liste_moyenne_decalage_de_l_origine (n,dj):   # Calcule la liste des valeurs l'écart quadratique moyen pour n valeurs de j séparées de dj
    L_moy_dec_origine=[]
    for j in range (dj*c,dj*c*n,dj):
        L_moy_dec_origine.append(moyenne_decalage_de_l_origine (j))
    return (L_moy_dec_origine)

X = np.arange(dj*c*dt,dj*c*dt*n,dj*dt)
L_moy_dec_origine = liste_moyenne_decalage_de_l_origine (n,dj)


plt.plot (X,L_moy_dec_origine,'b+',label = "Gamma={}e{} ; Potentiel=sinusoïdal ; Np={}e{} ; T={}K " .format(v_gamma,e_gamma,v_pas,e_pas,T))

plt.xlabel("Temps (s)",fontsize=15)
plt.ylabel("Ecart quadratique moyen (m^2)",fontsize=15)
plt.title("Ecart quadratique moyen en fonction du temps", fontsize=18)
plt.legend(prop={'size': 10}, title="Paramètres de simulation", loc="upper right")
plt.grid()

plt.show()
    






# X=[]
# for i in range(N):
#     print(i)
#     print(t)
#     t=t+dj
#     x=0
#     for j in range(0,pas-t):
#         for k in range(P):
#             x+=(atome[k].x[j+t]-atome[k].x[j])**2
#     X.append(x/P/(pas-t))

# lr = scipy.stats.linregress(T,X)
# print(lr[0]/2,dt,np.log(gamma)/np.log(10))

# Xtrue=T*lr[0]
# plt.scatter(T,X,label="P= "+str(P))
# plt.scatter(T,Xtrue)

# plt.show()