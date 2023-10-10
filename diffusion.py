from functions import *
import random as rd
import numpy as np
import matplotlib.pyplot as plt
import scipy

V=Vsin
Vprime=Vsinprime
P=1
N=20
dj=10000
t=0
T=np.arange(t,t+N*dj,dj)*dt
print(dt,np.log(gamma)/np.log(10))
simulation_RPMD(P,V,Vprime)



X=[]
for i in range(N):
    print(i)
    print(t)
    t=t+dj
    x=0
    for j in range(0,pas-t):
        for k in range(P):
            x+=(atome[k].x[j+t]-atome[k].x[j])**2
    X.append(x/P/(pas-t))

lr = scipy.stats.linregress(T,X)
print(lr[0]/2,dt,np.log(gamma)/np.log(10))

Xtrue=T*lr[0]
plt.scatter(T,X,label="P= "+str(P))
plt.scatter(T,Xtrue)

plt.show()