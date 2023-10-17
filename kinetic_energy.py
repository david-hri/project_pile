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


L_Ec = list_kinetic_energy(N)

M_Ec=int(np.mean(L_Ec)*10)/10
std=np.std(L_Ec)

IC=int(1.96*std/np.sqrt(N)*10)/10

std=int(std*10)/10
    
X1 = np.arange(N)

plt.plot (X1,L_Ec,'b+',label = "\u03B3={}e{} ; Potentiel=sinusoïdal ; Np={}e{} ; T={}K " .format(v_gamma,e_gamma,v_pas,e_pas,T))
plt.plot (X1,[M_Ec for i in range (N) ],'r',label="Energie cinétique moyenne={}K".format(M_Ec))
plt.plot (X1,[M_Ec-IC for i in range (N) ],'g',alpha=0.8,linewidth=1.2,label= "\u03C3={} ; IC={}".format(std,IC) )
plt.plot (X1,[M_Ec+IC for i in range (N) ],'g',alpha=0.8,linewidth=1.2)

plt.xlabel("Répétition",fontsize=15)
plt.ylabel("Energie cinétique moyenne (en K)",fontsize=15)
plt.title("Calcul de l'énergie cinétique moyenne", fontsize=18)
plt.legend(prop={'size': 10}, title="Paramètres de simulation et résultats", loc="upper right")
plt.grid()

plt.show()