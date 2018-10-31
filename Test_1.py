import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import matplotlib.cm as cm
T=200
B = 1.0/(T*constants.Boltzmann)
#print('B',B)
size=5
def E(s,J,size):
    Eng=0
    E_tot=0

    for i in range(size):
        for j in range(size):
            Eng+=s[i][j]*(s[(i+1)%size][j]+s[i][(j+1)%size]+s[(i-1)%size][j]+s[i][(j-1)%size])
    Eng=.5*J*Eng #.5 to remove double counting
    return -Eng
#def Norm(s,B):

numbers=range(-1,0) + range(1,2)
s= np.random.choice(numbers,(size,size))
CS =plt.matshow(s,vmin=-1, vmax=1,cmap="cool")
grid=True
#print (E(s,1,size))
s_bef=np.zeros_like(s)
s_aft=np.zeros_like(s)
s_bef=s
for q in range (100):
    for i in range (size):
        for j in range(size):
            E_bef = E(s_bef,1,size)
            s_aft[i][j]=-s[i][j]
            E_aft=E(s_aft,1,size)
            Del_E=E_aft-E_bef
            if (Del_E<0):
                '''p=np.exp(-B*Del_E)
                decide = np.random.ranf()
                if(decide<p):
                    s[i][j]=s_bef[i][j]
                else:'''
                s[i][j]=s_bef[i][j]
            
            elif(Del_E>0):
                s[i][j]=s_aft[i][j]
           
        
'''for i in range(size):
    for j in range (size):
        print (s[i][j], s_aft[i][j], s_bef[i][j])
      
        E_bef=E(s_bef,1,size)    
        if (s[i][j]==1):
            s_aft[i][j]=-1
            E_aft=E(s_aft,1,size)
            Del_E=E_aft-E_bef
        elif (s[i][j]==-1):
            s_aft[i][j]=1
            E_aft=E(s_aft,1,size)
            Del_E=E_aft-E_bef
        print Del_E
        if (Del_E<0):
            
            p=np.exp(-B*Del_E)
            decide = np.random.ranf()
            if(decide<p):
                s[i][j]=s_bef[i][j]
            else:
            s[i][j]=s_bef[i][j]
            print ("Cat")
        elif(Del_E>0):
            s[i][j]=s_aft[i][j]
            print ("Dog")
        print (s[i][j], s_aft[i][j], s_bef[i][j])
       
'''
    
CS =plt.matshow(s,vmin=-1, vmax=1,cmap="cool")
grid=True
#print (E(s,1,size))
plt.show()
