import sys
sys.modules[__name__].__dict__.clear()
import numpy as np
from numpy.linalg import eig, inv;
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d
import xlwt
import pandas as pd
import itertools as it


N = 10 # take multiple of 2 
# Fermi distribution function
def fermi(E,mu, T): #kb = 1boltzman constant
    if E > mu:
        return (math.exp(-(E - mu)/T))/(math.exp(-(E - mu)/T) + 1)
    else:
        return 1/(math.exp((E - mu)/T) + 1)
def savedata(Tm,Um,Deltam,name):
    book = xlwt.Workbook()
    sh = book.add_sheet('sheet')
    sh.write(0,0,'T')
    sh.write(0,1,'U')
    sh.write(0,2,'$\Delta$')
    
    for i in range(len(Tm)):
        sh.write(i+1,0,Tm[i])
        sh.write(i+1,1,Um[i])
        sh.write(i+1,2,Deltam[i])
    book.save(name)

#retrieing data from an excel file
def retrieve(name):
    file = pd.read_excel(name)
    T = file['T'];
    U = file['U'];
    Delta = file['$\Delta$'];
    Tm = [];Um = [];Deltam = []
    for i in range(len(T)):
            if  not math.isnan(float(T[i])):
                Tm.append(round(float(T[i]),3));
                Um.append(round(float(U[i]),5));
                Deltam.append(round(float(Delta[i]),3))
    return Tm,Um,Deltam

class Hubbard_Model:
    pass



#Finding The Energy Eigen Values
#basis c up d, c down,  
def Hamilt(tab, tbc, tca, U, mu, delta):
    H = [];
    for i in range(6*N*N):
        H.append([]);
        for j in range(6*N*N):
             if ((j == i + 2 or j == i - 4) and i in range(0,6*N*N,6) ) or ((j == i - 2 or j == i + 4) and i in range(2,6*N*N,6)):
                 H[i].append(-tab)
             elif ((j == i + 2 or j == i - 4) and i in range(1,6*N*N,6)) or ((j == i - 2 or j == i + 4) and i in range(3,6*N*N, 6)):
                 H[i].append(tab)
             elif (j == i + 4 or j == i - 6*N - 2) and i in range(0, 6*N*N,6) and ((i//(6*N))%2 == 0 ): 
                 H[i].append(-tca)
             elif  (j == i - 4 or j == i + 6*N - 4) and i in range(4,6*N*N,6) and ((i//(6*N))%2 == 0 ):
                 H[i].append(-tca)
             elif (j == i + 4 or j == i - 6*N + 4) and i in range(0, 6*N*N,6) and ((i//(6*N))%2 == 1 ):
                 H[i].append(-tca)
             elif (j == i - 4 or j == i + 6*N + 2) and i in range(4, 6*N*N,6) and ((i//(6*N))%2 == 1 ) :
                 H[i].append(-tca)
             elif (j == i + 4 or j == i - 6*N - 2) and i in range(1, 6*N*N,6) and ((i//(6*N))%2 == 0 ): 
                 H[i].append(tca)
             elif  (j == i - 4 or j == i + 6*N - 4) and i in range(5,6*N*N,6) and ((i//(6*N))%2 == 0 ):
                 H[i].append(tca)
             elif (j == i + 4 or j == i - 6*N + 4) and i in range(1, 6*N*N,6) and ((i//(6*N))%2 == 1 ):
                 H[i].append(tca)
             elif (j == i - 4 or j == i + 6*N + 2) and i in range(5, 6*N*N,6) and ((i//(6*N))%2 == 1 ) :
                 H[i].append(tca)
             elif  (j == i + 2 or j == i - 6*N + 2) and i in range(2, 6*N*N, 6) and ((i//(6*N))%2 == 0 ):  
                 H[i].append(-tbc)
             elif  (j == i + 2 or j == i - 6*N + 8) and i in range(2, 6*N*N, 6) and ((i//(6*N))%2 == 1 ):  
                H[i].append(-tbc)
             elif  (j == i - 2 or j == i + 6*N - 8) and i in range(4, 6*N*N, 6) and ((i//(6*N))%2 == 0 ):  
                H[i].append(-tbc)
             elif  (j == i - 2 or j == i + 6*N - 2) and i in range(4, 6*N*N, 6) and ((i//(6*N))%2 == 1 ):  
               H[i].append(-tbc)
             elif  (j == i + 2 or j == i - 6*N + 2) and i in range(3, 6*N*N, 6) and ((i//(6*N))%2 == 0 ):  
                 H[i].append(tbc)
             elif  (j == i + 2 or j == i - 6*N + 8) and i in range(3, 6*N*N, 6) and ((i//(6*N))%2 == 1 ):  
                H[i].append(tbc)
             elif  (j == i - 2 or j == i + 6*N - 8) and i in range(5, 6*N*N, 6) and ((i//(6*N))%2 == 0 ):  
                H[i].append(tbc)
             elif  (j == i - 2 or j == i + 6*N - 2) and i in range(5, 6*N*N, 6) and ((i//(6*N))%2 == 1 ):  
               H[i].append(tbc)
             elif j == i + 1 and i in range(0,6*N*N,2):
                 H[i].append(-U*delta)
             elif j == i - 1 and i in range(1,6*N*N,2):
                 H[i].append(-U*delta.conjugate())
             elif j == i and i in range(0,6*N*N,2):
                 H[i].append(-mu + U*(delta*delta.conjugate())) #-mu
             elif j == i and i in range(1,6*N*N,2):
                 H[i].append(mu + U*(delta*delta.conjugate())) #mu
             else:
                 H[i].append(0)
    for i in range(0,N,1): 
          if i%2 == 0 and i != 0 :
              #coorections for tca
              H[6*i*N][(i - 1)*6*N - 2] = 0
              H[6*i*N+1][(i - 1)*6*N - 1] = 0
              H[6*i*N-2][6*(i+1)*N] = 0
              H[6*i*N-1][6*(i+1)*N+1] = 0
              
              H[i*6*N][i*6*N - 2] = -tca #PBC
              H[i*6*N - 2][i*6*N] = -tca
              H[i*6*N + 1][i*6*N - 1] = tca
              H[i*6*N - 1][i*6*N+1] = tca
          if i%2 == 0:
              #H[6*(i+1)*N+4][6*(i+2)*N - 4] = 0 #redundant
              #H[6*(i+1)*N+5][6*(i+2)*N - 3] = 0 #redundant
              H[6*(i+2)*N - 3][6*(i+1)*N+5]= 0
              H[6*(i+2)*N - 4][6*(i+1)*N+4] = 0
              
              H[6*i*N+4][6*(i+1)*N-4] = 0
              #H[6*(i+1)*N-4][6*i*N+4] = 0 #redundant
              H[6*i*N+5][6*(i+1)*N-3] = 0
              #H[6*(i+1)*N-3][6*i*N+5] = 0 #redundant
              H[6*i*N+4][6*(i+2)*N-4] = -tbc
              H[6*(i+2)*N-4][6*i*N+4] = -tbc
              H[6*i*N+5][6*(i+2)*N-3] = tbc
              H[6*(i+2)*N-3][6*i*N+5] = tbc
         
          #PBC
          H[i*6*N][(i+1)*6*N - 4] = -tab
          H[(i+1)*6*N - 4][i*6*N] = -tab
          H[i*6*N+1][(i+1)*6*N-3] = tab
          H[(i+1)*6*N-3][i*6*N+1] = tab
          if i != 0:
             #corrections for tab
             H[i*6*N][i*6*N-4] = 0
             H[i*6*N-4][i*6*N] = 0
             H[i*6*N+1][i*6*N-3] = 0
             H[i*6*N-3][i*6*N+1] = 0 
              
             H[6*i][(N-1)*6*N + -2 + 6*i] = -tca
             H[(N-1)*6*N - 2 + 6*i][6*i] = -tca
             H[6*i+1][(N-1)*6*N + -1 + 6*i] = tca
             H[(N-1)*6*N - 1 + 6*i][6*i+1] = tca
          
          H[2+6*i][6*N*(N-1)+4 + 6*i] = -tbc
          H[6*N*(N-1)+4 +6*i][2+6*i] = -tbc
          H[3+6*i][6*N*(N-1)+5] = tbc
          H[6*N*(N-1)+5][3+6*i] = tbc
    H[0][6*N*N-2] = -tca
    H[6*N*N-2][0] = -tca
    H[1][6*N*N-1] = tca
    H[6*N*N-1][1] = tca
    return H


  
#band structure
def FBZ1():
    # half of FBZ to the full BZ
    Kx = []; Ky = [];
    for i in range(N):
        index = int(N*i - i*(i-1)/2);
        for j in range(N - i):
            Kx.append(-math.pi + 2*math.pi/(N-1)*i);
            Ky.append(-Kx[index]/math.sqrt(3) + 2*math.pi/math.sqrt(3) + j/(N-i)*(2*Kx[index]/math.sqrt(3) - 2*math.pi/math.sqrt(3)))
    for i in range(len(Kx)):
        Kx.append(-Kx[i]);
        Ky.append(-Ky[i]);
    return Kx, Ky
def FBZ2():
    Kx = []; Ky = []
    #one fourth of zone to the full
    for i in range(int(N/math.sqrt(6*math.sqrt(3)))):
        for j in range(int(N/math.sqrt(2*math.sqrt(3)))):
            Kx.append(2*math.pi/3/int(N/math.sqrt(6*math.sqrt(3)))  + (2*math.pi/3- 2*math.pi/3/( int(N/math.sqrt(6*math.sqrt(3))))  )*i/( int(N/math.sqrt(6*math.sqrt(3))) ));
            Ky.append(2*math.pi/math.sqrt(3) - 2*math.pi/math.sqrt(3)*j/( int(N/math.sqrt(2*math.sqrt(3))) ));
    for i in range(int(N/math.sqrt(6*math.sqrt(3)))):
        for j in range(int(N/math.sqrt(2*math.sqrt(3))) - round(math.sqrt(3)*i)):
            Kx.append(2*math.pi/3 + 2*math.pi/3*i/( int(N/math.sqrt(6*math.sqrt(3))) ));
            Ky.append((2*math.pi/math.sqrt(3) -i*2*math.pi/math.sqrt(3)/int(N/math.sqrt(6*math.sqrt(3))))*(1 - j/( int(N/math.sqrt(2*math.sqrt(3)))- round(math.sqrt(3)*i) )));
    for i in range(len(Kx)):
        Kx.append(Kx[i]);
        Ky.append(-Ky[i]);
    for i in range(len(Kx)):
        Kx.append(-Kx[i]);
        Ky.append(Ky[i])
    return Kx, Ky
def quasibandstructure(tab, tbc, tca,U, mu, delta, kx, ky):    
    # 6 x 6 matrix in k-space to find 
    HK = []
    for i in range(6):
        HK.append([])
        for j in range(6):
            HK[i].append(0)
    HK[0][0] = -mu + U*delta*delta.conjugate();
    HK[2][2] = HK[0][0]
    HK[4][4] = HK[0][0]
    
    HK[1][1] =  mu + U*delta*delta.conjugate();
    HK[3][3] =  HK[1][1]
    HK[5][5] =  HK[1][1]
    
    HK[0][1] =  U*delta
    HK[2][3] = HK[0][1];
    HK[4][5] = HK[0][1]
    
    HK[1][0] = U*delta.conjugate()
    HK[3][2] = HK[1][0]
    HK[5][4] = HK[1][0]
    
    HK[0][2] = -tab*complex(1+math.cos(kx), -math.sin(kx))  
    HK[2][0] = HK[0][2].conjugate()
    HK[1][3] =  tab*complex(1+math.cos(kx),  math.sin(kx)) 
    HK[3][1] = HK[1][3].conjugate();

    HK[0][4] = -tca*complex(1+math.cos(math.sqrt(3)/2*ky-kx/2),math.sin(math.sqrt(3)/2*ky-kx/2)) 
    HK[4][0] = HK[0][4].conjugate() 
    HK[1][5] =  tca*complex(1+math.cos(math.sqrt(3)/2*ky-kx/2),-math.sin(math.sqrt(3)/2*ky-kx/2))  
    HK[5][1] = HK[1][5].conjugate()
    
    HK[2][4] = -tbc*complex(1+math.cos(math.sqrt(3)/2*ky+kx/2),math.sin(math.sqrt(3)/2*ky+kx/2))
    HK[4][2] = HK[2][4].conjugate();
    HK[3][5] =  tbc*complex(1+math.cos(math.sqrt(3)/2*ky+kx/2),-math.sin(math.sqrt(3)/2*ky+kx/2))
    HK[5][3] = HK[3][5].conjugate()
    EigenValues, EigenVec = eig(HK) 
    EigenValues = EigenValues.real
    return list(EigenValues);

def bandstructure(tab, tbc, tca, mu, kx, ky):    
    # 6 x 6 matrix in k-space to find 
    
    HK = []
    for i in range(3):
        HK.append([])
        for j in range(3):
            HK[i].append(0)
    HK[0][0] = -mu 
    HK[1][1] = HK[0][0]
    HK[2][2] = HK[0][0]
    
    HK[0][1] = -tab*complex(1+math.cos(kx), -math.sin(kx))  
    HK[1][0] = HK[0][1].conjugate()

    HK[0][2] = -tca*complex(1+math.cos(math.sqrt(3)/2*ky-kx/2),math.sin(math.sqrt(3)/2*ky-kx/2)) 
    HK[2][0] = HK[0][2].conjugate() 
    
    HK[1][2] = -tbc*complex(1+math.cos(math.sqrt(3)/2*ky+kx/2),math.sin(math.sqrt(3)/2*ky+kx/2))
    HK[2][1] = HK[1][2].conjugate();
    EigenValues, EigenVec = eig(HK) 
    EigenValues = EigenValues.real
    return list(EigenValues);
   
tab = 1
tbc = 1
tca = 1
mu = 0
U = 3
Energy = []
candlistdelta = [i/1000 for i in range(0,200,10)] #set of candidate delta values 
delta = []    
#band structure
x, y = FBZ2();
Kx3 = [];Ky3 = [];
Kx6 = [];Ky6 = [];
for i in range(len(x)):
    Kx3.append(x[i]);Kx3.append(x[i]);Kx3.append(x[i]);
    Kx6.append(x[i]);Kx6.append(x[i]);Kx6.append(x[i]); Kx6.append(x[i]);Kx6.append(x[i]);Kx6.append(x[i])
    Ky3.append(y[i]);Ky3.append(y[i]);Ky3.append(y[i])
    Ky6.append(y[i]);Ky6.append(y[i]);Ky6.append(y[i]); Ky6.append(y[i]);Ky6.append(y[i]);Ky6.append(y[i])
quasibandvalues = []
bandvalues = []
for i in range(len(x)):
    #values1 = quasibandstructure(tab, tbc, tca,U,  mu, 0 ,  x[i], y[i]) #put delta value 
    values2 = bandstructure(tab, tbc, tca,  mu ,  x[i], y[i]) 
    #quasibandvalues = quasibandvalues + values1 # calling standard data
    bandvalues = bandvalues + values2
#selecting band values from EigenValues    
def pickupband(tab, tbc, tca, U, mu, delta, EigenValues, n = 6, Kx = Kx6, Ky = Ky6):
   index = []
   newEV = list(EigenValues);
   band = []
   if n == 3:
       band = bandvalues
   elif n == 6:
       band = quasibandvalues
   elif n == 0:
       for i in range(len(Kx)):
           values = bandstructure(tab, tbc, tca,  mu ,  Kx[i], Ky[i]) 
           band = band + values
   for i in range(len(band)):
       index.append(0);
       minn = 10**3
       for k in range(6*N*N):
           if abs(band[i] - EigenValues[k]) < minn and EigenValues[k] in newEV:
              index[i] = k
              minn = abs(band[i] - EigenValues[k])
       newEV.remove(EigenValues[index[i]]) 
       if len(newEV) == 0:
           break;
   return list(EigenValues[index])



# printing textbook type (2D) band structure diagram (delta = 0 perturbation off):
#H = Hamilt(tab, tbc, tca, U, mu, 0)
#EigenValues, EigenVec = eig(H)
#EigenValues = EigenValues.real    
##
##
#Kx = []; step = []
#Ky = []; 
#for i in range(int(N/math.sqrt(3))):
#    Kx.append(0 + 2*math.pi/3*i/int(N/math.sqrt(3)))   
#    Ky.append(4*math.pi/math.sqrt(3) - i*2*math.pi/math.sqrt(3)/int(N/math.sqrt(3)))
#for i in range(int(N/2/math.sqrt(3)) ):
#    Kx.append(2*math.pi/3 - i*2*math.pi/3/int(N/2/math.sqrt(3)));
#    Ky.append(2*math.pi/math.sqrt(3))
#EnergyValues = pickupband(tab, tbc, tca, U, mu, 0 , EigenValues, n = 0,Kx = Kx,Ky = Ky)
#for i in range(len(Kx)):
#    step.append(i);step.append(i);step.append(i);
#plt.scatter(step, EnergyValues, color = 'blue')




#Printing 3D band dispersion


#H = Hamilt(tab, tbc, tca, U, mu, 0)
#EigenValues, EigenVec = eig(H)
#EigenValues = EigenValues.real  # calling our matrix data
#ObsValues = pickupband(tab, tbc, tca, U, mu, 0 ,EigenValues, n = 3,Kx = Kx3,Ky = Ky3)
#fig = plt.figure(figsize = (10, 7))
#ax = plt.axes(projection ="3d")
#ax.scatter(Kx3, Ky3,bandvalues, color = 'red')
#ax.scatter(Kx3,Ky3,ObsValues, color = 'blue')
#plt.legend(['Standard curve'],['EigenValues of the matrix'])
#ax.set_xlabel('kx')
#ax.set_ylabel('ky')
#ax.set_zlabel('Energy')
#plt.title('N = 25')
#plt.show()    





#Code to produce phase diagram of the model and print in excel sheet
T = 0.7
Um = [i/100 for i in range(200,300,10)]
for m in range(len(Um)):
   Energy.append([])
   for k in range(len(candlistdelta)):
      H = Hamilt(tab, tbc, tca, Um[m], mu, candlistdelta[k])
      EigenValues, EigenVec = eig(H)
      EigenValues = EigenValues.real
      #at 0 temperature, H - muN has zero expectation value in the ground state
      # minizing energy method
      temp = 0
      for i in range(6*N*N):
           temp = temp + EigenValues[i]*fermi(EigenValues[i],0,T);
      Energy[m].append(temp)
   delta.append(candlistdelta[Energy[m].index(min(Energy[m]))])
Tm = [T]*len(Um);
a, b, c = retrieve('C:/Users/sarve/OneDrive/Desktop/New folder/PHY665/Hubbard.xls')
Tm =  a + Tm;
Um = b + Um;
delta = c + delta;
savedata(Tm,Um,delta,'C:/Users/sarve/OneDrive/Desktop/New folder/PHY665/Hubbard.xls')


#Code to find Delta versus T
#T = [i/100 for i in range(50, 51, 1)]; #Tc is between 0.9 and 1.1. for all ts = 1, mu=0, Tc = 0.95, U = 3
#for m in range(len(T)):
#   Energy.append([])
#   for k in range(len(candlistdelta)):
#      H = Hamilt(tab, tbc, tca, U, mu, candlistdelta[k])
#      EigenValues, EigenVec = eig(H)
#      EigenValues = EigenValues.real
#      #at 0 temperature, H - muN has zero expectation value in the ground state
#      # minizing energy method
#      temp = 0
#      for i in range(6*N*N):
#           temp = temp + EigenValues[i]*fermi(EigenValues[i],0,T[m]);
#      Energy[m].append(temp)
#   delta.append(candlistdelta[Energy[m].index(min(Energy[m]))])



