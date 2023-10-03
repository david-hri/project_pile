import random as rd
import numpy as np
import matplotlib.pyplot as plt

class Atome:

    def __init__(self):
            self.x = [0,0,0]
            self.v = [0,0,0]



m=1.6726*10**(-27)
ev=1.6*10**(-19)
a=0.3*10**(-10)
V0=4.8e-20
kb=1.38*10**(-23)
hb=6.626*10**(-34)/2
hb/=np.pi
e=a/100
w0 = np.sqrt((8*V0/(a**2*m)))
P=1
T=2000
pas=300000
wP=P*kb*T/hb
dt=10**(-17)
gamma=10e15

i=0
K=P**2*m*kb**2*T**2/hb**2/100
atome={}

def Vp(x):
    return V0*(1-(x/a)**2)**2

def Vpprime (x) :
    V_prime = (4*V0/a**2) * ((x**3/a**2) - x)
    return (V_prime)

def Vsin(x):
    return V0*np.sin(2*np.pi*x/a)**2

def Vsinprime(x):
    return 4*V0*np.pi/a*np.sin(2*np.pi*x/a)*np.cos(2*np.pi*x/a)



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

def main_DM(Vprime):
    F = Vprime(atome[0].x[-1]) - m*gamma*atome[0].v[-1] + np.sqrt ((2*m*gamma*kb*T)/dt)*rd.gauss(0,1)
    pos = 2*atome[0].x[-1] - atome[0].x[-2]+dt**2 * F/m
    vit = (3*atome[0].x[-1]-4*atome[0].x[-2]+atome[0].x[-3])/(2*dt)
    atome[0].x.append(pos)
    atome[0].v.append(vit)
    
def simulation_DM(Vprime): #main function which computes the position and the speed of one particle in a given distribution of potential
    init(1)
    for i in range(pas):
        if i%100000==0:
            print(i)
        main_DM(Vprime)

def main_RPMD(etape,N,Vprime):  #integration of the motion equation, 
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


