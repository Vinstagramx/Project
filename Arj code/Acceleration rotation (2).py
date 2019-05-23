
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Q = np.genfromtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/nonInertialdata.csv', delimiter = ",", skip_header = 1)

accelraw=Q[:,[1,2,3]] #raw accel data from csv

gyroraw=Q[:,[7,8,9]] #raw gyro data from csv

gravraw=Q[:,[4,5,6]]

def pitch(gx):
    M=[[1,0,0],[0,sp.cos(gx),-sp.sin(gx)],[0,sp.sin(gx),sp.cos(gx)]]
    return M #rotations in the pitch, roll and yaw diretions respectively in matrix form
    
def roll(gy):
    M=[[sp.cos(gy),0,sp.sin(gy)],[0,1,0],[-sp.sin(gy),0,sp.cos(gy)]]
    return M

def yaw(gz):
    M=[[sp.cos(gz),-sp.sin(gz),0],[sp.sin(gz),sp.cos(gz),0],[0,0,1]]
    return M

def totalrot(x,y,z): # just a function I defined to multiply three matrices together cos I couldnt find one on scipy
    A=np.matmul(x,y)
    B=np.matmul(A,z)
    return B 

antipitchset = []
antirollset = []
antiyawset=[]

for x in range (0,len(gyroraw[:,0])):
    M=pitch(-gyroraw[x,0])
    antipitchset.append(M) 
    M=roll(-gyroraw[x,1])
    antirollset.append(M)    #creates a set of matrices which will invert the pitch on each of the readings
    M=yaw(-gyroraw[x,2])
    antiyawset.append(M)
                            # so when I apply the first matrix in the list to the first gyro vector it will invert it to be with respect to ground


    
totrotset=[] # this becomes the list of all transformation matrices which will return the accel to the ground frame in all three axes (pitch roll and yaw)
 
for x in range(0,len(gyroraw[:,0])):
    M=totalrot(antipitchset[x],antirollset[x],antiyawset[x])
    totrotset.append(M)

acclist=[] # rotated accel into ground frame

for x in range(0,len(gyroraw[:,0])):
    M=np.matmul(totrotset[x],accelraw[x])
    acclist.append(M)

accelreal=np.asarray(acclist) #turned the rotated data into an array to plot it

gravlist=[]
for x in range(0,len(gyroraw[:,0])):
    M=np.matmul(totrotset[x],gravraw[x])
    gravlist.append(M)

gravreal=np.asarray(gravlist)

#plt.plot(sp.arange(0,len(accelreal[:,0])),accelreal[:,0],)
#plt.title('accel in x')
#plt.show()
#plt.plot(sp.arange(0,len(accelreal[:,1])),accelreal[:,1],) #dirty accel plots after rotation
#plt.title('accel in y')
#plt.show()
#plt.plot(sp.arange(0,len(accelreal[:,2])),accelreal[:,2],)
#plt.title('accel in z')
#plt.show()



#%%
rang=0.01*sp.arange(0,11550) # the time data on x axis of graph in seconds
from scipy.signal import lfilter # some random filter package I found on internet
n = 50  # the larger n is, the smoother curve will be
b = [1.0 / n] * n
a = 1


xx = sp.signal.savgol_filter(accelreal[:,0],501,3)
#plt.plot(rang, xx, linewidth=2, linestyle="-") # plots of xyz clean accels
#plt.title('clean accel in x direction')
#plt.show()
yy = sp.signal.savgol_filter(accelreal[:,1],501,3)
#plt.plot(rang,yy,linewidth=2,linestyle="-")
#plt.title('clean accel in y direction')
#plt.show()
zz = sp.signal.savgol_filter(accelreal[:,2],501,3)
#plt.plot(rang,zz,linewidth=2,linestyle="-")
#plt.title('clean accel in z direction')
#plt.show()
#%%
t=Q[:,0][5240:8800]-[Q[0,0]]*3560
#az=accelreal[:,2]
#
#delta_t=[]
#for x in range(0,t.size -1):
#    m=t[x+1]-t[x]
#    delta_t.append(m)
#
#n=t.size
#dt=np.mean(delta_t)
#freq=np.fft.fftfreq(n,dt)  # a nice helper function to get the frequencies  
#coeff=np.fft.fft(az) # returns the coefficients for each frequency
#
#df=10 # accel changes must happen in a timespan more than 0.1s to be physical
#env=np.exp(-(freq/(2*df))**2)
#filt=env*coeff
#
#plt.plot(freq,coeff)
#plt.plot(freq,filt)
#plt.show()
#Az=np.fft.ifft(filt)
#plt.plot(t,Az)

#%%

def fftfilter(t,data):
    del_t=[] 
    for x in range(0,len(t)-1):
        m=t[x+1]-t[x]
        del_t.append(m)
        
    n=len(t)
    dt=np.mean(del_t) # the data from MEMs has slightly varied dt 
    freq=np.fft.fftfreq(n,dt)  # gives fourier frequencies  
    coeff=np.fft.fft(data) #fourier coefficients
    
    df=2 # accel changes must happen in a timespan more than 0.25s to be physical
    env=np.exp(-(freq/(2*df))**2) # gaussian of std 10 in fourier space to filter high freq
    filt=env*coeff # filtered fourier coefficients
    
    data_filt=np.fft.ifft(filt)#filtered data
    return [freq,coeff,filt,data_filt] #returns the useful data
    













#%%
#from scipy.integrate import cumtrapz
#
#Az2=np.real(Az[4500:9000])
#Az2_inv=Az2[::-1]
#vz=cumtrapz(Az2)
#vz_inv=cumtrapz(Az2_inv)
#
#wv=[]
#for x in range (0,vz.size-1):
#    w=(vz.size-x)/vz.size
#    m=vz[x]*w
#    wv.append(m)
#wv=np.asarray(wv)
#
#wv_inv=[]
#for x in range (0,vz_inv.size-1):
#    w=(vz_inv.size-x)/vz_inv.size
#    m=vz_inv[x]*w
#    wv_inv.append(m)
#wv_inv=np.asarray(wv_inv)
#
#wv_noninv=wv_inv[::-1]
#
#WV=[]
#for x in range(0,wv.size-1):
#    k=np.mean([wv[x],wv_noninv[x]])
#    WV.append(k)
#    
#    
#    
    
    
    
#%%
from scipy.integrate import cumtrapz

def boundary(X):
    #inverting the list and integrating fwds and bkwds
    Xinv=X[::-1]
    V=cumtrapz(X,initial=0)
    Vinv=cumtrapz(Xinv,initial=0)
    #creating a list of the weighted fwds and bkwds integrated values
    wV=[]
    for n in range (0,V.size):
        w=((V.size-1-n)/V.size-1)
        m=V[n]*w
        wV.append(m)
    wV=np.asarray(wV)
    
    wVinv=[]
    for n in range (0,Vinv.size):
        w=((Vinv.size-1-n)/Vinv.size-1)
        m=Vinv[n]*w
        wVinv.append(m)
    wVinv=np.asarray(wVinv)
    
    wVnoninv=wVinv[::-1] #swapping inverted list right way round again
    wVmean=[]
    for x in range(0,V.size):
        k=np.mean([wV[x],wVnoninv[x]])
        wVmean.append(k)
    wVmean=np.asarray(wVmean)
    return wVmean



#%%
ax=accelreal[:,0][4500:9000]
ay=accelreal[:,1][4500:9000]
az=accelreal[:,2][4500:9000]
t_trim=t[4500:9000]

ax2=fftfilter(t_trim,ax)[3]
ay2=fftfilter(t_trim,ay)[3]
az2=fftfilter(t_trim,az)[3]



SX=boundary(boundary(ax2))
SY=boundary(boundary(ay2))
SZ=boundary(boundary(az2))


from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(SX,SY,SZ)

#%%
