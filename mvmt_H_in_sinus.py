from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt

simulation_RPMD(1,Vsin,Vsinprime)  

def tracé_potentiel(Xp):
    V=[]
    for x in Xp :
        V.append(Vsin(x))
    return V

X1 = [i*dt for i in range(len(atome[0].x))]
Xp = np.arange(-4*a,4*a,8*a/pas)
V = tracé_potentiel(Xp)
   
fig, ax1 = plt.subplots()
    
ax1.plot(atome[0].x,X1)
ax1.set_xlabel("Position (m)")
ax1.set_ylabel("Temps (s)")
#fig.title("Position du proton dans un potentiel sinusoidale en fonction du temps")
    
ax2=ax1.twinx()
ax2.set_ylabel("potentiel sinusoidale (J)")
ax2.plot(Xp,V, 'r')
ax2.set_ylim(0,5e-20)
    
plt.title ("Evolution d'un proton dans un potentiel sinusoidal")
plt.show()