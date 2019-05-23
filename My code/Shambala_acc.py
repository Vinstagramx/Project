import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import numpydoc as npd
import sphinx as sph
import nose as ns
import pykalman as pk

Data = np.loadtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/nonInertialdata.csv', delimiter = ",", skiprows = 1)

rawacc = Data[:, [1,2,3]]
rawgyro = Data[:, [7,8,9]]

def pitch(a):
    M = [[1,0,0],[0,sp.cos(a),-sp.sin(a)],[0,sp.sin(a),sp.cos(a)]]
    return M #rotations in the pitch, roll and yaw diretions respectively in matrix form
    
def roll(b):
    N = [[sp.cos(b),0,sp.sin(b)],[0,1,0],[-sp.sin(b),0,sp.cos(b)]]
    return N

def yaw(c):
    O = [[sp.cos(c),-sp.sin(c),0],[sp.sin(c),sp.cos(c),0],[0,0,1]]
    return O

def multiply_three_matrices(x,y,z): # just a function I defined to multiply three matrices together cos I couldnt find one on scipy
    Q = np.matmul(np.matmul(x,y),z)
    return Q


invertedpitch = []
invertedroll = []
invertedyaw = []
rotation = []
acceleration_ground_list = []
for x in range (0,len(rawgyro[:,0])):
    M = pitch(-rawgyro[x,0])
    invertedpitch.append(M)  
    #creates a set of matrices which will invert the pitch on each of the readings
    N = roll(-rawgyro[x,1])
    invertedroll.append(N)  
                        # so when I apply the first matrix in the list to the first gyro vector it will invert it to be with respect to ground
    O = roll(-rawgyro[x,2])
    invertedyaw.append(O) 
    
    Q = multiply_three_matrices(invertedpitch[x],invertedroll[x],invertedyaw[x])
    rotation.append(M)
    
    R=np.matmul(rotation[x],rawacc[x])
    acceleration_ground_list.append(R)                               # this becomes the list of all transformation matrices which will return the accel to the ground frame in all three axes (pitch roll and yaw)

 # rotated accel into ground frame
    

acceleration_ground = np.asarray(acceleration_ground_list)

plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,0])
plt.title("Acceleration in the x-direction relative to the ground")
plt.show()

plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,1]) #dirty accel plots after rotation
plt.title('Acceleration in the y-direction relative to the ground')
plt.show()
plt.plot(sp.arange(0,len(rawgyro[:,0])), acceleration_ground[:,2])
plt.title('Acceleration in the z-direction relative to the ground')
plt.show()
