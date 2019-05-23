import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

gps_data=np.genfromtxt('C:/Users/Arjun Kolli/Documents/Summer Project/Raw Data/Stealth/gpsData.csv', delimiter = ",", skip_header = 1)
lat=gps_data[:,1] -[gps_data[0,1]]*len(gps_data[:,1])
long=gps_data[:,2] - [gps_data[0,2]]*len(gps_data[:,2])
Z=gps_data[:,3] 

baro_data=np.genfromtxt('C:/Users/Arjun Kolli/Documents/Summer Project/Raw Data/Stealth/barometicAltitudeData.csv', delimiter = ",", skip_header = 1)
P=1000*baro_data[:,1]
h=baro_data[:,2]

#%%

#def height_from_pressure2(P):
#    Tb=288.15
#    Lb=-0.0065
#    Pb=101325
#    g=9.80665
#    R=8.3144598
#    M=0.0289644
#    
#    K=(P/Pb)**(-(R*Lb)/(g*M))
#    h=(Tb/Lb)*(K-1)
#    return h

def height_from_pressure(P):
    Tb=288.15
    Lb=-0.0065
    Pb=101325
    g=9.80665
    R=8.3144598
    M=0.0289644
    K=sp.log(P/Pb)
    C=(-R*Tb)/(g*M)
    h=C*K
    return h
    
h2=[]
for x in range (0,120): 
    H=height_from_pressure(P[x])
    h2.append(H)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot(sx2[0:3498:159],sy2[0:3498:159],h[62:84:1])
#plt.show()

fig = plt.figure()
Ax=fig.add_subplot(111,projection='3d')
Ax.plot(lat,long,h[0:116])
plt.show()

#%%

