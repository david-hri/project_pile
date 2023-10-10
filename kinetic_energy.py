import matplotlib.pyplot as plt
import numpy as np
from functions import *

N = 10

def list_kinetic_energy(N):
    L_Ec = []
    for i in range (N) :
        simulation_RPMD(P,Vsinprime)
        L_Ec.append(E_c())
        print(i+1)
    return (L_Ec)

def av_kinetic_energy (L_Ec):
    Av=0
    for i in range(N):
        Av += L_Ec[i]
    Av_Ec = Av/N
    return (Av_Ec)


L_Ec = list_kinetic_energy(N)

Av_Ec=int(av_kinetic_energy (L_Ec)*10e2)/10e2
    
X1 = np.arange(N)

plt.plot (X1,L_Ec,'b+',label = "Gamma=1e14 ; Potentiel=sinusoïdal ; Np=5e6 ; T={}K " .format(T))
plt.plot (X1,[Av_Ec for i in range (N) ],'r',label="Energie cinétique moyenne={}".format(Av_Ec))

plt.xlabel("Répétition",fontsize=15)
plt.ylabel("Energie cinétique moyenne (en K)",fontsize=15)
plt.title("Calcul de l'énergie cinétique moyenne", fontsize=18)
plt.legend(prop={'size': 10}, title="Paramètres de simulation", loc="upper right")
plt.grid()

plt.show()
        
    