
from functions import *

simulation_RPMD(P,Vp,Vpprime)
plt.plot([i for i in range(len(atome[0].x))],atome[0].x)
print(atome[0].x[-1])
plt.show()


