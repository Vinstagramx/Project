import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import numpydoc as npd
import sphinx as sph
import nose as ns
import pykalman as pk
from scipy.integrate import cumtrapz

Data = np.loadtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/nonInertialdata.csv', delimiter = ",", skiprows = 1)
g = 9.80665
rawacc = Data[:, [1,2,3]] * g
rawgyro = Data[:, [7,8,9]] * g

def rot180(theta): #rotates matrix by 180 degrees
    return np.matrix([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])

def rotation(x, y, z):
    '''
    rotation matrix with Tait-Bryan angles
    where x, y and z are rotation angles about X, Y and Z
    rotation performed in order Z -> X -> Y
    https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles
    '''
    c1 = np.cos(z)
    s1 = np.sin(z)
    c2 = np.cos(x)
    s2 = np.sin(x)
    c3 = np.cos(y)
    s3 = np.sin(y)
    matrix = np.array([[((c1 * c3) - (s1 * s2 * s3)), -c2 * s1, ((c1 * s3) + (c3 * s1 * s2))], [((c3 * s1) + (c1 * s2 * s3)), c1 * c2, ((s1 * s3) - (c1 * c3 * s2))], [-c2 * s3, s2, c2 * c3]])
    rot = rot180(np.deg2rad(180)) * matrix
    rotatedacc = np.matmul(rawacc,rot)
    return rotatedacc

acceleration_ground_list = []
for x in range (0,len(rawgyro[:,0])):
    acceleration_ground_list.extend(rotation(rawgyro[x,0], rawgyro[x,1], rawgyro[x,2]))

# this becomes the list of all transformation matrices which will return the accel to the ground frame in all three axes (pitch roll and yaw)
  
# rotated accel into ground frame
acceleration_ground = np.asarray(acceleration_ground_list)

plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,0])
plt.title("Acceleration in the x-direction relative to the ground/ms^-2")
plt.savefig(xacc)
plt.show()
plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,1])
plt.title("Acceleration in the y-direction relative to the ground/ms^-2")
plt.savefig(yacc)
plt.show()
plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,2])
plt.title("Acceleration in the z-direction relative to the ground/ms^-2")
plt.savefig(zacc)
plt.show()
