import numpy as np
import matplotlib.pyplot as plt
import random as rd
from scipy.optimize import curve_fit

### Valeur numérique des paramètres :

## Paramètres physique :

# Masse du proton
m_p = 1.67e-27  

# Constante de Boltzman 
k_b = 1.38e-23 # J K-1

# Température
T = 2000 # K


# Forme du potentiel 
V_0 = 4.8e-20 # J
a = 3e-11 # m
 
def potentiel (x) :
    V = V_0 * (np.sin((2*np.pi*x)/a))**2
    return V

def dérivée_potentiel (x) :
    V_prime = (4*V_0*np.pi/a) *np.sin((2*np.pi*x)/a) * np.cos((2*np.pi*x)/a)
    return (V_prime)
    
# Tracé du potentiel 

def tracé_potentiel(Xp):
    V=[]
    for x in Xp :
        V.append(potentiel(x))
    return V
 
# Fréquence d'ocsillation du proton dans le potentiel 
omega_O = np.sqrt((8*V_0/(a**2*m_p)))


## Paramètres de simulation :

# Temps caractéristique de simulation 
delta_t = 1e-17

# Gamma
gamma = 4e10


# Nombre de points
nbr_p = 8000000


### Algorithme de résolution

def algorithme_verlet (n):
    
    # Initialisation 
    L_pos = [0]                     
    L_vit = [0]
    F = - dérivée_potentiel (L_pos[0]) - m_p*gamma*L_vit[0] + np.sqrt ((2*m_p*gamma*k_b*T)/delta_t)*rd.gauss(0,1)
    
    # Calcul de la position, de la vitesse et de la force à l'instant delta_t
    pos = L_pos[0] + delta_t * L_vit[0] + delta_t**2 * F/(2*m_p)
    L_pos.append(pos)
    
    vit = (3*L_pos[1]-4*L_pos[0]+L_pos[0])/(2*delta_t)  
    L_vit.append(vit)
    
    F = - dérivée_potentiel (L_pos[1]) - m_p*gamma*L_vit[1] + np.sqrt ((2*m_p*gamma*k_b*T)/delta_t)*rd.gauss(0,1)
    
    
    for i in range (1,n-1) :      
        
        # Calcul de la position à l'instant t+delta_t
        pos = 2*L_pos[i] - L_pos[i-1]+delta_t**2 * F/m_p           
        L_pos.append(pos)
        
        # Calcul de la vitesse à l'instant t+delta_t
        vit = (3*L_pos[i+1]-4*L_pos[i]+L_pos[i-1])/(2*delta_t)        
        L_vit.append(vit)
        
        # Calcul de la force à l'instant t+delta_t
        F = - dérivée_potentiel (L_pos[i+1]) - m_p*gamma*L_vit[i+1] + np.sqrt ((2*m_p*gamma*k_b*T)/delta_t) *rd.gauss(0,1)
        
    return L_pos,L_vit

## Calcul de l'energie cinétique moyenne

def energie_cinétique_moyenne_rec (L_vit, nbr_p) :           # calcul de l'intégrale avec la méthode des rectangles
    V = 0
    for i in range (len(L_vit)) :
        V += L_vit[i]**2                                    
    E_c = m_p/(2*nbr_p) *V
    return E_c

## Méthode des décalages de l'origine

def moyenne_decalage_de_l_origine (L_pos,j) :   # Calcule la valeur moyenne du décalage de l'origine pour j fixé
    S = 0
    for k in range (nbr_p-j):
        S+=(L_pos[j+k]-L_pos[k])**2
    Moy_dec = S/(nbr_p-j)
    return Moy_dec

def  liste_moyenne_decalage_de_l_origine (n,delta_j,L_pos):   # Calcule la liste des valeurs moyennes du décalage de l'origine pour n valeurs de j séparées de delta_j
    L_moy_dec_origine=[]
    for j in range (delta_j*n*8,delta_j*n*9,delta_j):
        L_moy_dec_origine.append(moyenne_decalage_de_l_origine (L_pos,j))
    return (L_moy_dec_origine)

def tirage_liste_moyenne_decalage_de_l_origine(k):          # Cette foncttion permet de calculer une liste de listes contenant l'écart quadratique moyen obtenu grâce à la méthode du décalage de l'origine pour n valeurs de j séparées de delta_j. Chaque Liste est réalisé sur k positions différentes calculés grâce a l'algorithme de Verlet
    tir_L_moy_dec_origine = []
    for i in range (k):
        L_pos=algorithme_verlet (nbr_p)[0]
        tir_L_moy_dec_origine.append(liste_moyenne_decalage_de_l_origine (n,delta_j,L_pos))
    return tir_L_moy_dec_origine

def moyenne_liste_moyenne_decalage_de_l_origine(k):
    moy_L_moy_dec_origine = []
    tir_L_moy_dec_origine = tirage_liste_moyenne_decalage_de_l_origine(k)
    S = 0
    M = 0
    for j in range (0,n) :
        for i in range (k):
            S += tir_L_moy_dec_origine[i][j]
        M = S/k
        moy_L_moy_dec_origine.append(M) 
        S = 0
    return (moy_L_moy_dec_origine)


# Calcul avec le modèle Scipy

# def func (x,a,b):
#     return a*x+b
        
        
# def droite_diffusion(param,X2):       #Fonction qui permet de tracer la droite de l'écart quadratique moyen en fonction du temps ayant pour pente le coefficient de diffusion moyen
#     CD,O = param[0],param[1]
#     print(param)
#     print(CD,O)
#     Ec_quad_moy = []
#     for j in X2 :
#         Ec_quad_moy.append ((2*CD*j)+O)
#     return Ec_quad_moy

def regression_lineaire(X, Y):
    X = np.array(X)
    Y = np.array(Y)
    A = np.vstack([X, np.ones(len(X))]).T
    m, c = np.linalg.lstsq(A, Y, rcond=None)[0]
    return m, c

def droite_diffusion(m, c, X):
    return m * X + c

def coefficient_de_diffusion (moy_tir_L_moy_dec_origine, P0,X):                 #Calcul le coefficient de diffusion partir de la l-ième liste de l'écart quadratique moyen en fonction du temps grâce au module scipy
    m,c = regression_lineaire(X,moy_tir_L_moy_dec_origine)
    return (m*1/2,c)

def capcaplameilleure (p):
    L=[]
    while p<=10000:
        L.append("capcap est trop forte")
        p+=1
    return (L)

if __name__=='__main__' :
    
    print ("ATTENTION VALEUR GAMMA LEGENDE GRAPHIQUE")
    
    # Lsave=np.load("save.npy")
    # L_moy_dec_origine=list(Lsave)
    
    n=20
    delta_j = 1000
    
    k=10
    
    L_pos, L_vit = algorithme_verlet(nbr_p)
    
    X = np.arange(delta_j*delta_t*n*5+delta_j*delta_t,delta_j*delta_t*n*6,delta_j*delta_t)

    moy_tir_L_moy_dec_origine = moyenne_liste_moyenne_decalage_de_l_origine(k)
    
    print (len(moy_tir_L_moy_dec_origine),len(X))
    D, c = regression_lineaire(X, moy_tir_L_moy_dec_origine)
    Y = droite_diffusion(D, c, X)
    
    CD  = int(D*10**9)/1e9
    print(CD)

    plt.plot(X,moy_tir_L_moy_dec_origine,'b+')
    plt.plot(X,Y,'r',label = "Gamma=4e10 ; CD={} ; k={}" .format(CD,k))

    plt.xlabel("Temps (s)",fontsize=15)
    plt.ylabel("Ecart quadratique moyen (m^2)",fontsize=15)
    plt.title("Ecart quadratique moyen en fonction du temps", fontsize=18)
    plt.grid()
    plt.legend(prop={'size': 12})
    plt.show()
    
    # print (coefficient_de_diffusion (tir_L_moy_dec_origine,0, P0))
    # print (coefficient_de_diffusion (tir_L_moy_dec_origine,1, P0))
    # print (coefficient_de_diffusion (tir_L_moy_dec_origine,2, P0))
    # print (coefficient_de_diffusion (tir_L_moy_dec_origine,3, P0))
    # print(moyenne_coefficient_de_diffusion(tir_L_moy_dec_origine,k))

    
    # Lsave=np.array(L_moy_dec_origine)
    # print(len(Lsave))
    # np.save('save',Lsave)


#print (capcaplameilleure(1))
