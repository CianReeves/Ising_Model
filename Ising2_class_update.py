import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from numpy import random
from scipy.spatial import Delaunay
random.seed(69)
N=50 #Number of different temperature points taken for certain results
T=1.0#Initial temperature
iterations=300 #number of iterations after the lattice is warmed up
Temp=np.zeros(N) 
Mag_avg=np.zeros(N) #List to store average magnetization for certain temps
Mag_avg_sqr=np.zeros(N) #Stores second moment of magnetization
suscept=np.zeros(N) #stores susceptibility
E_avg,E_avg_sqr,spec_heat=np.zeros(N),np.zeros(N),np.zeros(N)
class Ising:
    def __init__(self,size,iterations,h,T):
        self.size = size
        self.s=np.random.choice((-1,1),(size,size))
        self.iterations=iterations #Number of iterations to warm up lattice
        self.J=1
        self.K=constants.Boltzmann
        self.T=T*self.J/self.K
        self.B=1.0/(self.T*self.K)
        self.h=h
    
    def Rand_lattice(self):
        
        grid = random.rand(30,2)
        tri=Delaunay(grid)
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(grid[:,0],grid[:,1],'bo', ms=2)
        ax.triplot(grid[:,0],grid[:,1],tri.simplices)
        ax.set_aspect('equal')
        pindex=3
        print tri.vertex_neighbor_vertices[1]
        print grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][1][0]
       
        print grid[3][1]
        
        #if tri.neighbors>=0:
        #print grid[tri.simplices]
        #print grid[tri.simplices][1][2]
        #print grid
    
    def Neighbors_rand(self):
        num_points=(self.size)**2
        grid=np.array([[a+(b%2)*.5,b*np.sqrt(3)/2] for a in np.linspace(0,self.size-1,self.size) for b in np.linspace(0,self.size-1,self.size)])
        #grid = random.rand(num_points,2)
        tri=Delaunay(grid)
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(grid[:,0],grid[:,1],'bo')
        ax.triplot(grid[:,0],grid[:,1],tri.simplices)
        ax.set_aspect('equal')
        print tri.neighbors
        neigh=[]
        coord=[]
        #if tri.neighbors>=0:
        #pindex = 5
        for pindex in range(num_points):
            temp_dist=0
            nearest=[]#list which holds the nearest neighbours of the lattice point pindex
            
            j=0
            for i in range (len(tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]])):#loops through all the points of the lattice
                x_dist= grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][0]-grid[pindex][0]
                y_dist= grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][1]-grid[pindex][1]
                
                dist = float("{0:.2f}".format(np.sqrt(x_dist**2+y_dist**2)))
                if (dist<temp_dist or temp_dist==0):#checks if the point is the current nearest point 
                    nearest.append(grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][0])
                    nearest.append(grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][1])
                    repeats=[]
                if (dist==temp_dist):#this checks if there is multiple points that are all the nearest neighbours 
                    repeats.append(grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][0])
                    repeats.append(grid[tri.vertex_neighbor_vertices[1][tri.vertex_neighbor_vertices[0][pindex]:tri.vertex_neighbor_vertices[0][pindex+1]]][i][1])
                temp_dist=dist
            nearest.append(repeats) #combines nearest initial nearest neighbour with any repeated nearest neighbours
            neigh.append(nearest) #combines all nearest neighbours into one list
            coord.append(grid[pindex])
        print neigh
        print coord
        #print grid[tri.simplices][1][1]
      #  print grid
        
        
        
        
        
    def Triangular(self):
        grid=np.array([[a+(b%2)*.5,b*np.sqrt(3)/2] for a in np.linspace(0,self.size-1,self.size) for b in np.linspace(0,self.size-1,self.size)])
        tri=Delaunay(grid)
        
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(grid[:,0],grid[:,1],'bo')
        ax.triplot(grid[:,0],grid[:,1],tri.simplices)
        ax.set_aspect('equal')
        print tri.points
        
    def Neighbours_sqr(self,s,i,j): #Nearest Neighbours for square lattice
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
        for i in range(self.iterations): #Calculates magnetization for a series of warmed up lattices and then averages them
            Mag_avg[j]+=self.Mag(s_warm)
            Mag_avg_sqr[j]+=(self.Mag(s_warm))**2
            s_warm=self.Flip(s_warm)
        Mag_avg[j] = abs(Mag_avg[j])/self.iterations
        Mag_avg_sqr[j]=Mag_avg_sqr[j]/self.iterations
        suscept[j]=(1.0/self.T)*(Mag_avg_sqr[j]-Mag_avg[j]**2)
        return suscept
        
    def Energy(self,s_warm):
        E_eq=0
        for i in range(self.size):
            for j in range(self.size):
                E_eq+=s_warm[i][j]*self.Neighbours(s_warm,i,j)
        E_eq=E_eq/(2*(self.size)**2)#factor of 2 to remove double counting of energies.
        return -E_eq
    
    def Specific_heat(self,s_warm,j):
        for i in range(self.iterations):
            E_avg[j]+=self.Energy(s_warm)
            E_avg_sqr[j]+=self.Energy(s_warm)**2
            s_warm=self.Flip(s_warm)
        E_avg[j] = abs(E_avg[j])/self.iterations
        E_avg_sqr[j]=E_avg_sqr[j]/self.iterations
        spec_heat[j]=(self.J**2/(self.K**2*self.T**2))*(E_avg_sqr[j]-E_avg[j]**2)
        return spec_heat
        
        
        
        
    


'''f1=plt.figure()
    ax=f1.add_subplot(131)
    im=plt.imshow(latt.s,vmin="-1",vmax="1",cmap="plasma")'''   
'''for j in range(N):
    latt=Ising(10,300,0,T) #Initializes Ising class with different temperatures
    for i in range(100): #Warms up Lattice
        s_warm=latt.Flip(latt.s)
    #susceptibility=latt.susceptibility(s_warm,j)
    Specific_heat=latt.Specific_heat(s_warm,j)
    Temp[j]=T
    T+=.05'''
latt=Ising(300,300,0,2)

latt.Neighbors_rand()
'''f3=plt.figure() #figure for average magnetization.
ax3=f3.add_subplot(111)
ax3.plot(Temp,Mag_avg,'ro')

f4=plt.figure()#Figure for susceptibility
ax4=f4.add_subplot(111)
ax4.plot(Temp,susceptibility,'b*')'''


'''f5=plt.figure()
ax5=f5.add_subplot(111)
ax5.plot(Temp,Specific_heat,'g*')

f6=plt.figure()
ax6=f6.add_subplot(111)
ax6.plot(Temp,E_avg,'b*')'''
plt.show()
