import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from numpy import random
#random.seed(69)
N=50 #Number of different temperature points taken for certain results
T=1.0#Initial temperature
iterations=300 #number of iterations after the lattice is warmed up
Temp=np.zeros(N) 
Mag_avg=np.zeros(N) #List to store average magnetization for certain temps
Mag_avg_sqr=np.zeros(N) #Stores second moment of magnetization
suscept=np.zeros(N) #stores susceptibility
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
        
    
    def Neighbours(self,s,i,j): #Nearest Neighbours for square lattice
        R=s[(i+1)%self.size][j]
        L=s[(i-1)%self.size][j]
        U=s[i][(j+1)%self.size]
        D=s[i][(j-1)%self.size]
        return(R+L+U+D)
    
    def dE(self,i,j): #Calclates del_E
        Eng=self.s[i][j]*self.Neighbours(self.s,i,j)
        return -2.0*Eng*self.J +2*self.h*self.s[i][j]
    
    def Flip(self,s): #Runs through lattice and flips each particles spin and decides whether or not to accept
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
   
    def Mag(self,s_warm): #Calclates the magnetization
        Mag=0
        
        for i in range (self.size):
            for j in range(self.size):
                Mag+=s_warm[i][j]
        
        Mag=1.0/(self.size**2)*Mag
        
        return Mag

    def susceptibility(self,s_warm,j): #Calculates susceptibility at a certain temperature
            for i in range(iterations): #Calculates magnetization for a series of warmed up lattices and then averages them
                Mag_avg[j]+=self.Mag(s_warm)
                Mag_avg_sqr[j]+=(self.Mag(s_warm))**2
                s_warm=self.Flip(s_warm)
            Mag_avg[j] = abs(Mag_avg[j])/iterations
            Mag_avg_sqr[j]=Mag_avg_sqr[j]/iterations
            suscept[j]=(1.0/self.T)*(Mag_avg_sqr[j]-Mag_avg[j]**2)
            return suscept
        
''' 
        
    def Energy(self,s_warm):
        E_eq=0
        for i in range(size):
            for j in range(size):
                E_eq+=s_warm[i][j]*self.Neighbours(s_warm,i,j)
        E_eq=E_eq/2'''
    


'''f1=plt.figure()
    ax=f1.add_subplot(131)
    im=plt.imshow(latt.s,vmin="-1",vmax="1",cmap="plasma")'''   
for j in range(N):
    latt=Ising(10,100,0,T) #Initializes Ising class with different temperatures
    for i in range(100): #Warms up Lattice
        s_warm=latt.Flip(latt.s)
    susceptibility=latt.susceptibility(s_warm,j)
    Temp[j]=T
    T+=.05
    
f3=plt.figure() #figure for average magnetization.
ax3=f3.add_subplot(111)
ax3.plot(Temp,Mag_avg,'ro')

f4=plt.figure()#Figure for susceptibility
ax4=f4.add_subplot(111)
ax4.plot(Temp,susceptibility,'b*')
plt.show()
