import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from numpy import random
random.seed(69)

class Ising:
    def __init__(self,size,warm_up,h,T):
        self.size = size
        self.s=np.random.choice((-1,1),(size,size))
        self.warm_up=warm_up #Number of iterations to warm up lattice
        self.J=1
        self.K=constants.Boltzmann
        self.T=T*self.J/self.K
        self.B=1.0/(self.T*self.K)
        self.h=h
    
    def Neighbours(self,s,i,j):
        R=s[(i+1)%self.size][j]
        L=s[(i-1)%self.size][j]
        U=s[i][(j+1)%self.size]
        D=s[i][(j-1)%self.size]
        return(R+L+U+D)
    
    def dE(self,i,j):
        Eng=self.s[i][j]*self.Neighbours(self.s,i,j)
        return -2.0*Eng*self.J +2*self.h*self.s[i][j]
    
    def Flip(self,s):
        for i in range (self.size):
            for j in range(self.size):
                s[i][j]=-s[i][j]#
                Del_E=self.dE(i,j)
                if (Del_E>0.0):
                    p=np.exp(-self.B*Del_E)
                    if(np.random.ranf()<p):
                        s[i][j]=s[i][j]
                    else:
                        s[i][j]=-s[i][j]
                
                elif(Del_E<=0.0):
                    s[i][j]=s[i][j]
        return s
   
    def Mag(self,s_warm):
        Mag=0
        
        for i in range (self.size):
            for j in range(self.size):
                Mag+=s_warm[i][j]
        
        Mag=1.0/(self.size**2)*Mag
        
        return Mag

    def Energy(self,s_warm):
        E_eq=0
        for i in range(size):
            for j in range(size):
                E_eq+=s_warm[i][j]*self.Neighbours(s_warm,i,j)
        E_eq=E_eq/2
N=50
T=1.0
Mag_avg=np.zeros(N)
Temp=np.zeros(N)

for j in range(N):
    latt=Ising(10,100,0,T)


    '''f1=plt.figure()
    ax=f1.add_subplot(131)
    im=plt.imshow(latt.s,vmin="-1",vmax="1",cmap="plasma")'''               
    
    for i in range(100): #Warms up Lattice
        s_warm=latt.Flip(latt.s)

    #ax=f1.add_subplot(132)
    #im=plt.imshow(s_warm,vmin="-1",vmax="1",cmap="plasma")

    for i in range(200): #Calculates magnetization for a series of warmed up lattices and then averages them
        Mag_avg[j]+=latt.Mag(s_warm)
        s_warm=latt.Flip(s_warm)
    Mag_avg[j] = abs(Mag_avg[j])/100 
    Temp[j]=T
    T+=.05

f=plt.figure()
ax=f.add_subplot(111)
ax.plot(Temp,Mag_avg,'ro')
plt.show()
