import Sphinx 
import numpydoc 
import nose
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pykalman

Q = np.genfromtxt('C:/Users/Arjun Kolli/Documents/Summer Project/Raw Data/Stealth/stealthraw.csv', delimiter = ",", skip_header = 1)

accelraw=Q[:,[1,2,3]] #raw accel data from csv

gyroraw=Q[:,[7,8,9]] #raw gyro data from csv

def pitch(gx):
    M=[[1,0,0],[0,sp.cos(gx),-sp.sin(gx)],[0,sp.sin(gx),sp.cos(gx)]]
    return M #rotations in the pitch, roll and yaw diretions respectively in matrix form
    
def roll(gy):
    M=[[sp.cos(gy),0,sp.sin(gy)],[0,1,0],[-sp.sin(gy),0,sp.cos(gy)]]
    return M

def yaw(gz):
    M=[[sp.cos(gz),-sp.sin(gz),0],[sp.sin(gz),sp.cos(gz),0],[0,0,1]]
    
def totalrot(x,y,z): # just a function I defined to multiply three matrices together cos I couldnt find one on scipy
    A=np.matmul(x,y)
    B=np.matmul(A,z)
    return B 

antipitchset = []

for x in range (0,len(gyroraw[:,0])):
    M=pitch(-gyroraw[x,0])
    antipitchset.append(M)      #creates a set of matrices which will invert the pitch on each of the readings
                                # so when I apply the first matrix in the list to the first gyro vector it will invert it to be with respect to ground
antirollset = []

for x in range (0,len(gyroraw[:,1])):
    M=roll(-gyroraw[x,1])
    antirollset.append(M)

antiyawset=[]

for x in range (0,len(gyroraw[:,2])):
    M=roll(-gyroraw[x,2])
    antiyawset.append(M) 
    
totrotset=[] # this becomes the list of all transformation matrices which will return the accel to the ground frame in all three axes (pitch roll and yaw)
 
for x in range(0,len(gyroraw[:,0])):
    M=totalrot(antipitchset[x],antirollset[x],antiyawset[x])
    totrotset.append(M)

acclist=[] # rotated accel into ground frame

for x in range(0,len(gyroraw[:,0])):
    M=np.matmul(totrotset[x],accelraw[x])
    acclist.append(M)

accelreal=np.asarray(acclist) #turned the rotated data into an array to plot it


plt.plot(sp.arange(0,len(accelreal[:,0])),accelreal[:,0],)
plt.title('accel in x')
plt.show()
plt.plot(sp.arange(0,len(accelreal[:,1])),accelreal[:,1],) #dirty accel plots after rotation
plt.title('accel in y')
plt.show()
plt.plot(sp.arange(0,len(accelreal[:,2])),accelreal[:,2],)
plt.title('accel in z')
plt.show()



rang=0.01*sp.arange(0,11550) # the time data on x axis of graph in seconds
from scipy.signal import lfilter # some random filter package I found on internet
n = 50  # the larger n is, the smoother curve will be
b = [1.0 / n] * n
a = 1

xx = lfilter(b,a,accelreal[:,0])
plt.plot(rang, xx, linewidth=2, linestyle="-") # plots of xyz clean accels
plt.title('clean accel in x direction')
plt.show()
yy = lfilter(b,a,accelreal[:,1])
plt.plot(rang,yy,linewidth=2,linestyle="-")
plt.title('clean accel in y direction')
plt.show()
zz = lfilter(b,a,accelreal[:,2])
plt.plot(rang,zz,linewidth=2,linestyle="-")
plt.title('clean accel in z direction')
plt.show()















