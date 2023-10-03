from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt




simulation_RPMD(P,Vp,Vpprime)
txt="simu_P="+str(P)+"T="+str(T)+"pas="+str(pas)+".txt"
with open(txt, 'w') as fichier:
    for i in range(P):
        phrase=""
        print(i)
        for j in atome[i].x:
            phrase+=str(j)+" "
        
        fichier.write(phrase+"\n")
        phrase=""
        for j in atome[i].v:
            phrase+=str(j)+" "
        fichier.write(phrase+"\n")


 