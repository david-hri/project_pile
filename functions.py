import random as rd
import numpy as np
import matplotlib.pyplot as plt

class Atome:

    def __init__(self):
            self.x = [a,a,a]
            self.v = [0,0,0]



m=1.6726*10**(-27)
ev=1.6*10**(-19)
a=0.3*10**(-10)
V0=3*ev
kb=1.380649 *10**(-23)
hb=6.626*10**(-34)/2
hb/=np.pi
e=a/100
w0 = np.sqrt((8*V0/(a**2*m)))
P=1
T=2000
pas=10000000
wP=P*kb*T/hb
dt=1/1000*2*np.pi/w0
dt=10**(-17)
gamma=1/(1000*dt)

i=0
K=P**2*m*kb**2*T**2/hb**2/100
atome={}

def Vp(x):
    return V0*(1-(x/a)**2)**2

def Vpprime (x) :
    V_prime = (4*V0/a**2) * ((x**3/a**2) - x)
    return (V_prime)
def Vsin(x):
    return V0**2*np.sin(2*np.pi*x/a)**2

def Vsinprime(x):
    return 4*np.pi/a*V0**2*np.sin(2*np.pi*x/a)*np.cos(2*np.pi*x/a)



def init(N):   #creates N atoms with the same caracteristics
    for i in range(N):
        atome[i]=Atome()

def voisins(i,N): # gives the neighbours of a givin atom
    if N==1:
        return (0,0)
    if i==N-1:
        return (N-2,0)
    elif i==0:
        return (N-1,1)
    return i-1,i+1


def main_RPMD(etape,N,V,Vprime):  #integration of the motion equation, 
    for i in range(N):
        F = - 1/N*Vprime(atome[i].x[-1]) - m*gamma*atome[i].v[-1] + np.sqrt ((2*m*gamma*kb*T)/dt)*rd.gauss(0,1)
        pp,pn=voisins(i,N)
        Fk=-K*(2*atome[i].x[2+etape]-atome[pp].x[2+etape]-atome[pn].x[2+etape])
        F+=Fk
        pos = 2*atome[i].x[-1] - atome[i].x[-2]+dt**2 * F/m
        vit = (3*atome[i].x[-1]-4*atome[i].x[-2]+atome[i].x[-3])/(2*dt)
        atome[i].x.append(pos)
        atome[i].v.append(vit)

def simulation_RPMD(N,V,Vprime): #main function which computes the positions of the P particles in a given distribution of potential
    init(N)
    for i in range(pas):
        if i%100000==0:
            print(i)
        main_RPMD(i,N,V,Vprime)


def Veff(B,V):
    x=0
    for i in B:
        x+=V(i)/kb/T
    return x/len(B)

def Vmoy_for_P(N,V):
    Vmoy=[0]
    for i in range(pas):
        Vmoy.append((Vmoy[-1]*(i+1)+Veff([atome[j].x[-1] for j in range(N)],V))/(i+2))
    return Vmoy[-1]

#test coucou
