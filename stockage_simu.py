from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt




simulation_RPMD(P,Vpprime)
file="mvmt_particule n°{} ; Potentiel=sinusoïdal ; gamma={}e{} ; Np={}e{} ; T={}K ; P={}" .format(num_sim,v_gamma,e_gamma,v_pas,e_pas,T,P)

format_str = "%.4e"                                 # Cela permet ne sauvegarder qu'avec 4 chiffres significatifs
np.savetxt(file,atome[0].x,fmt=format_str)

L = np.loadtxt(file)
print(L)

# txt="simu_P="+str(P)+"T="+str(T)+"pas="+str(pas)+".txt"
# with open(txt, 'w') as fichier:
#     for i in range(P):
#         phrase=""
#         print(i)
#         for j in atome[i].x:
#             phrase+=str(j)+" "
        
#         fichier.write(phrase+"\n")
#         phrase=""
#         for j in atome[i].v:
#             phrase+=str(j)+" "
#         fichier.write(phrase+"\n")


 