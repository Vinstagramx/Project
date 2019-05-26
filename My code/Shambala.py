import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import numpydoc as npd
import sphinx as sph
import nose as ns
import pykalman as pk
import pandas as pd
from scipy.integrate import cumtrapz
from mpl_toolkits import mplot3d

'''
    Obtaining data and placing into arrays
'''
Data = np.loadtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/nonInertialdata.csv', delimiter = ",", skiprows = 1)
g = 9.80665
rawacc = Data[:, [1,2,3]] * g #multiply by g to correct acc
rawgrav = Data[:,[4,5,6]]
rawgyro = Data[:, [7,8,9]]
rawax = rawacc[:,0]
raway = rawacc[:,1]
rawaz = rawacc[:,2]
roll = rawgyro[:,0]
pitch = rawgyro[:,1]
yaw = rawgyro[:,2]
gx = rawgrav[:,0]
gy = rawgrav[:,1]
gz = rawgrav[:,2]
'''
    Rotating the acceleration to ground frame
'''
def rot180(theta): #rotates matrix by 180 degrees
    return np.matrix([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])

def taitbryan(x, y, z):
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
    matrix = np.matrix([[((c1 * c3) - (s1 * s2 * s3)), -c2 * s1, ((c1 * s3) + (c3 * s1 * s2))], [((c3 * s1) + (c1 * s2 * s3)), c1 * c2, ((s1 * s3) - (c1 * c3 * s2))], [-c2 * s3, s2, c2 * c3]])
    return rot180(np.deg2rad(180)) * matrix

specificrot = []
xacc = []
yacc = []
zacc = []
for x in range(0,len(gx)):

    rotation = taitbryan(pitch[x],roll[x],yaw[x])
    acc = rotation* np.matrix([[rawax[x]],[raway[x]],[rawaz[x]]])
    xacc.append(acc[0,0])
    yacc.append(-acc[1,0])
    zacc.append(acc[2,0])

# this becomes the list of all transformation matrices which will return the accel to the ground frame in all three axes (pitch roll and yaw)
# rotated accel into ground frame


'''
    Plotting unfiltered accelerations
'''
plt.plot(xacc)
plt.title("Acceleration in the x-direction relative to the ground/ms^-2")
plt.savefig('xacc.png', dpi = 500)
plt.show()
plt.plot(yacc)
plt.title("Acceleration in the y-direction relative to the ground/ms^-2")
plt.savefig('yacc.png', dpi = 500)
plt.show()
plt.plot(zacc)
plt.title("Acceleration in the z-direction relative to the ground/ms^-2")
plt.savefig('zacc.png', dpi = 500)
plt.show()

'''
    Trimming data
'''
time = Data[:,0][9500:25500]-[Data[0,0]]*16000

'''
    Fast Fourier Filter
'''
def fftfilter(time,data):
    delta_t=[] 
    for x in range(0,len(time)-1):
        m = time[x+1] - time[x]
        delta_t.append(m)
        
    n=len(time)
    dt=np.mean(delta_t) # the data from MEMs has slightly varied dt 
    freq=np.fft.fftfreq(n,dt)  # gives fourier frequencies  
    coeff=np.fft.fft(data) #fourier coefficients
    
    df=1 # accel changes must happen in a timespan more than 0.25s to be physical
    env=np.exp(-(freq/(2*df))**2) # gaussian of std 10 in fourier space to filter high freq
    filt=env*coeff # filtered fourier coefficients
    
    data_filt=np.fft.ifft(filt)#filtered data
    return [freq,coeff,filt,data_filt] #returns the useful data

"""
trimming and filtering the acceleration data
"""
filteredxacc = np.real(fftfilter(time,xacc[9500:25500])[3])
filteredyacc = np.real(fftfilter(time,yacc[9500:25500])[3])
filteredzacc = np.real(fftfilter(time,zacc[9500:25500])[3])

plt.plot(filteredxacc)
plt.title("Filtered acceleration in the x-direction relative to the ground/ms^-2")
plt.savefig('fxacc.png', dpi = 500)
plt.show()
plt.plot(filteredyacc)
plt.title("Filtered acceleration in the y-direction relative to the ground/ms^-2")
plt.savefig('fyacc.png', dpi = 500)
plt.show()
plt.plot(filteredzacc)
plt.title("Filtered acceleration in the z-direction relative to the ground/ms^-2")
plt.savefig('fzacc.png', dpi = 500)
plt.show()

"""
    Defining any acceleration below a certain value as being equal to zero - as rate of change is
    stil present in the random changes/noise in the very small accelerations - which affects the integration.
"""
for i, val in enumerate(filteredxacc):
    if val < 0.25:
        filteredxacc[i] = 0
        
for i, val in enumerate(filteredyacc):
    if val < 0.25:
        filteredyacc[i] = 0
        
for i, val in enumerate(filteredzacc):
    if val < 0.25:
        filteredzacc[i] = 0

"""
    Plotting velocities
"""

plt.plot(xvel)
plt.title("Velocity in the x-direction relative to the ground/ms^-1")
plt.savefig('fxvel.png', dpi = 500)
plt.show()
plt.plot(yvel)
plt.title("Velocity in the y-direction relative to the ground/ms^-1")
plt.savefig('fyvel.png', dpi = 500)
plt.show()
plt.plot(zvel)
plt.title("Velocity in the z-direction relative to the ground/ms^-1")
plt.savefig('fzvel.png', dpi = 500)
plt.show()


'''
    Integration function with boundary condition
'''
def intbound(X,t):
    #inverting the list and integrating fwds and bkwds
    X_flip = X[::-1]
    vel = cumtrapz(X,t,initial=0)
    vel_flip = cumtrapz(X_flip,t,initial=0)
    #creating a list of the weighted fwds and bkwds integrated values
    weighted_vel1 = []
    for n in range (0,len(vel)):
        w = ((len(vel)-1-n)/len(vel)-1) +1 #w=weight
        weighted_vel1.append(vel[n]*w)

    weighted_vel1 = np.asarray(weighted_vel1)
    
    weighted_vel_flip = []
    for n in range (0,len(vel_flip)):
        w = ((len(vel_flip)-1-n)/len(vel_flip)-1) +1
        weighted_vel_flip.append(vel_flip[n]*w)

    weighted_vel_flip = np.asarray(weighted_vel_flip)
    
    weighted_vel2 = weighted_vel_flip[::-1] #swapping inverted list right way round again

    totalvel = []
    for x in range(0,len(vel)):
        totalvel.append(weighted_vel1[x] + weighted_vel2[x])
        
    totalvel = np.asarray(totalvel)
    return totalvel


"""
Integration - twice from acc to position and plotting
"""
xvel = intbound(filteredxacc,time)
xpos = intbound(xvel,time)
yvel = intbound(filteredyacc,time)
ypos = intbound(yvel,time)
zvel = intbound(filteredzacc,time)
zpos = intbound(zvel,time)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot3D(xpos, ypos, zpos, 'gray')

'''
    Obtaining barometer data
'''
barometer_data = np.loadtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/barometicAltitudeData.csv', delimiter = ",", skiprows = 1)
pressure = barometer_data[:,1]*1000
phoneest_height = barometer_data[:,2]

def pres2alt(pressure):
    '''
    Determines altitude (in metres above sea level) from site pressure (in Pascals).

    Assumptions:
    ============================   ================
    Parameter                      Value
    ============================   ================
    Base pressure                  101325 Pa
    Temperature at zero altitude   288.15 K
    Gravitational acceleration     9.80665 m/s^2
    Lapse rate                     -6.5E-3 K/m
    Gas constant for air           287.053 J/(kgK)
    Relative Humidity              0%
    ============================   ================
    '''
    alt = 44331.5 - 4946.62 * pressure ** (0.190263)
    return alt
altitude = pres2alt(pressure)
plt.plot(altitude,phoneest_height)
plt.title("Altitude from phone against calculated altitude")
plt.savefig('alt.jpg', dpi = 500)
plt.show()
"""
    We will use the calculated altitude 
    data as the phone estimated altitude 
    is relative to the initial point (and thus is not the true altitude)
"""

'''
    Pressure assimilation
'''
trimmedalt = np.asarray(altitude[26:270])

"""
trimming the barometer data to the correct time
"""

def veltopos(vel,pos,t): #should return correct position in 1D. use boundary first to get vel from acc

    vel_right = []    #trimmed lists of vel to right of position point for each position point
    vel_left = [] #left of position point
    time_right = []
    time_left = []
    n = round(len(vel)/len(pos)) #ratio between measurement freq. tells you how many vel to skip per position - has to be an integer
    m = round(len(time)/len(pos))
    for x in range (0,len(pos)):
        vel_right.extend([vel[(x*n):len(pos)]])
        vel_left.extend([vel[0:(x*n)]])
        time_right.extend([time[(x*m):len(pos)]])
        time_left.extend([time[0:(x*m)]])

    right = []#forward trajectory to right of position with x0=x_position for each a point
    right_flip = []#backward right
    left = []
    left_flip = []
    for x in range (0,len(pos)):
        right.append(cumtrapz(vel_right[x],time_right,initial=pos[x]))

        right_flip_backwards = cumtrapz(vel_right[x][::-1],time_right,initial=0)
        right_flip.append(right_flip_backwards[::-1])

        left.append(cumtrapz(vel_left[x],time_left,initial=0))

        left_flip_backwards = cumtrapz(vel_left[x][::-1],time_left,initial=pos[x])
        left_flip.append(left_flip_backwards[::-1])     #gotta flip the backwards integrals cos the integrated data will be given bkwds to the fwd data 
   

    left_weight = [] # weighted trajectory of fwd and back for ONE altitude point (to the right)
    left_weight_tot = [] #list of all weighted trajectories to left of position points 
    right_weight = []
    right_weight_tot = []
     
    for t in range(0,len(pos)):
        m = len(left[t]) #length of a specific left traj for a specific position point
        k = len(right[t]).size

        if m==1:
            left_weight.append(left[t][0][0])#do this to avoid dividing by zero in the weighting process; cos if m=1, m-1=0
        else:
            for x in range (0,m):
                left_f_w = ((m-x-1)/(m-1))*left[t][0][x] #fwd weighted so most confident at points near position point
                left_b_w = (x/(m-1))*left_flip[t][0][x] #bkwd weighted so most confident at end start of tot traj
                left_weight.append(np.mean([left_f_w,left_b_w])) #add each weighted point to the traj to make up one whole traj
        left_weight_tot.append(left_weight)#add specific weighted traj for specif position point to set of all left traj
        if k==1:
            right_weight.append(right[t][0][0])
        else:
            for x in range (0,k):
                right_f_w = ((k-x-1)/(k-1))*right[t][0][x]
                right_b_w = (x/(k-1))*right_flip[t][0][x]
                right_weight.append(np.mean([right_f_w,right_b_w]))  
                #print(np.mean([right_f_w,right_b_w]))
        right_weight_tot.append(right_weight)#same deal but for traj to the right of position points 
        
    weight_tot = []#list of the total weighted trajectories for each position point
    for t in range (0,len(right_weight_tot)):
        weight_tot.append(left_weight_tot[t]+right_weight_tot[t])#just combine left and right trajectories
    full_trajectory = []#weight each trajectory, with highest confidence for the points near the position point
    for t in range(0,len(pos)):
        for x in range(0,len(weight_tot)):
            if x>=round(n*t): # this is weighting for points to right of position point
                w1 = (len(weight_tot)-1-x)/(len(weight_tot)-round(n*t)-1) #stuff to the left of position point
                full_trajectory.append(w1*weight_tot[t][x])                
            else:
                w2 = x/(round(n*t)-1)
                full_trajectory.append(w2*weight_tot[t][x])
                
    return full_trajectory
#%%
zpos=veltopos(zvel,trimmedalt,time)
plt.plot(zpos)

"""
    Converting latitude and longitude to x and y
"""
import math

__all__ = ['to_latlon', 'from_latlon']

K0 = 0.9996

E = 0.00669438
E2 = E * E
E3 = E2 * E
E_P2 = E / (1.0 - E)

SQRT_E = math.sqrt(1 - E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = _E * _E
_E3 = _E2 * _E
_E4 = _E3 * _E
_E5 = _E4 * _E

M1 = (1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256)
M2 = (3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
M3 = (15 * E2 / 256 + 45 * E3 / 1024)
M4 = (35 * E3 / 3072)

P2 = (3. / 2 * _E - 27. / 32 * _E3 + 269. / 512 * _E5)
P3 = (21. / 16 * _E2 - 55. / 32 * _E4)
P4 = (151. / 96 * _E3 - 417. / 128 * _E5)
P5 = (1097. / 512 * _E4)

R = 6378137

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"


def to_latlon(easting, northing, zone_number, zone_letter=None, northern=None, strict=True):
    """This function convert an UTM coordinate into Latitude and Longitude

        Parameters
        ----------
        easting: int
            Easting value of UTM coordinate

        northing: int
            Northing value of UTM coordinate

        zone number: int
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. More information see utmzones [1]_

        zone_letter: str
            Zone Letter can be represented as string values. Where UTM Zone
            Designators can be accessed in [1]_

        northern: bool
            You can set True or False to set this parameter. Default is None


       .. _[1]: http://www.jaworski.ca/utmzones.htm

    """
    if not zone_letter and northern is None:
        raise ValueError('either zone_letter or northern needs to be set')

    elif zone_letter and northern is not None:
        raise ValueError('set either zone_letter or northern, but not both')

    if strict:
        if not 100000 <= easting < 1000000:
            raise Exception('easting out of range (must be between 100.000 m and 999.999 m)')
        if not 0 <= northing <= 10000000:
            raise Exception('northing out of range (must be between 0 m and 10.000.000 m)')
    if not 1 <= zone_number <= 60:
        raise Exception('zone number out of range (must be between 1 and 60)')

    if zone_letter:
        zone_letter = zone_letter.upper()

        if not 'C' <= zone_letter <= 'X' or zone_letter in ['I', 'O']:
            raise Exception('zone letter out of range (must be between C and X)')

        northern = (zone_letter >= 'N')

    x = easting - 500000
    y = northing

    if not northern:
        y -= 10000000

    m = y / K0
    mu = m / (R * M1)

    p_rad = (mu +
             P2 * math.sin(2 * mu) +
             P3 * math.sin(4 * mu) +
             P4 * math.sin(6 * mu) +
             P5 * math.sin(8 * mu))

    p_sin = math.sin(p_rad)
    p_sin2 = p_sin * p_sin

    p_cos = math.cos(p_rad)

    p_tan = p_sin / p_cos
    p_tan2 = p_tan * p_tan
    p_tan4 = p_tan2 * p_tan2

    ep_sin = 1 - E * p_sin2
    ep_sin_sqrt = math.sqrt(1 - E * p_sin2)

    n = R / ep_sin_sqrt
    r = (1 - E) / ep_sin

    c = _E * p_cos**2
    c2 = c * c

    d = x / (n * K0)
    d2 = d * d
    d3 = d2 * d
    d4 = d3 * d
    d5 = d4 * d
    d6 = d5 * d

    latitude = (p_rad - (p_tan / r) *
                (d2 / 2 -
                 d4 / 24 * (5 + 3 * p_tan2 + 10 * c - 4 * c2 - 9 * E_P2)) +
                 d6 / 720 * (61 + 90 * p_tan2 + 298 * c + 45 * p_tan4 - 252 * E_P2 - 3 * c2))

    longitude = (d -
                 d3 / 6 * (1 + 2 * p_tan2 + c) +
                 d5 / 120 * (5 - 2 * c + 28 * p_tan2 - 3 * c2 + 8 * E_P2 + 24 * p_tan4)) / p_cos

    return (math.degrees(latitude),
            math.degrees(longitude) + zone_number_to_central_longitude(zone_number))


def latlon_xy(latitude, longitude, force_zone_number=None):
    """This function convert Latitude and Longitude to UTM coordinate

        Parameters
        ----------
        latitude: float
            Latitude between 80 deg S and 84 deg N, e.g. (-80.0 to 84.0)

        longitude: float
            Longitude between 180 deg W and 180 deg E, e.g. (-180.0 to 180.0).

        force_zone number: int
            Zone Number is represented with global map numbers of an UTM Zone
            Numbers Map. You may force conversion including one UTM Zone Number.
            More information see utmzones [1]_

       .. _[1]: http://www.jaworski.ca/utmzones.htm
    """
    if not -80.0 <= latitude <= 84.0:
        raise Exception('latitude out of range (must be between 80 deg S and 84 deg N)')
    if not -180.0 <= longitude <= 180.0:
        raise Exception('longitude out of range (must be between 180 deg W and 180 deg E)')

    lat_rad = math.radians(latitude)
    lat_sin = math.sin(lat_rad)
    lat_cos = math.cos(lat_rad)

    lat_tan = lat_sin / lat_cos
    lat_tan2 = lat_tan * lat_tan
    lat_tan4 = lat_tan2 * lat_tan2

    if force_zone_number is None:
        zone_number = latlon_to_zone_number(latitude, longitude)
    else:
        zone_number = force_zone_number

    lon_rad = math.radians(longitude)
    central_lon = zone_number_to_central_longitude(zone_number)
    central_lon_rad = math.radians(central_lon)

    n = R / math.sqrt(1 - E * lat_sin**2)
    c = E_P2 * lat_cos**2

    a = lat_cos * (lon_rad - central_lon_rad)
    a2 = a * a
    a3 = a2 * a
    a4 = a3 * a
    a5 = a4 * a
    a6 = a5 * a

    m = R * (M1 * lat_rad -
             M2 * math.sin(2 * lat_rad) +
             M3 * math.sin(4 * lat_rad) -
             M4 * math.sin(6 * lat_rad))

    easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000

    northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))

    if latitude < 0:
        northing += 10000000

    return easting, northing, zone_number


def latitude_to_zone_letter(latitude):
    if -80 <= latitude <= 84:
        return ZONE_LETTERS[int(latitude + 80) >> 3]
    else:
        return None


def latlon_to_zone_number(latitude, longitude):
    if 56 <= latitude < 64 and 3 <= longitude < 12:
        return 32

    if 72 <= latitude <= 84 and longitude >= 0:
        if longitude <= 9:
            return 31
        elif longitude <= 21:
            return 33
        elif longitude <= 33:
            return 35
        elif longitude <= 42:
            return 37

    return int((longitude + 180) / 6) + 1


def zone_number_to_central_longitude(zone_number):
    return (zone_number - 1) * 6 - 180 + 3

gps_data = np.loadtxt('C:/Users/user/Documents/Project/Roller Coaster Data/Shambala/gpsData.csv', delimiter = ",", skiprows = 1)
latitude = gps_data[:,1]
longitude = gps_data[:,2]

"""
    Trimming latitude and longitude data
"""
lattrim = latitude[93:245]
longtrim = latitude[93:245]

xy = []
origin = latlon_xy(lattrim[0],longtrim[0])
for x in range(0, len(lattrim) - 1):
    converted = latlon_xy(lattrim[x],longtrim[x])
    xy.append([converted[0] - origin[0], converted[1] - origin[1]])
    
xyTranspose = list(map(list, zip(*xy)))
x = xyTranspose[0]
y = xyTranspose[1]
"""
xyTranspose is a 1d list - [0] is latitude, [1] is longitude, [2] is zone number(useless)
access each element by [0][x] - gives xth element of latitude

"""

    


 