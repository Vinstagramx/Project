import numpy as np
import pandas as pd

AZ=np.real(fftfilter(t,accelreal[:,0][5240:8800])[3])#trimmed data we want to work with
H=np.asarray(h[55:92])
#
#t2=baro_data[:,0][55:95]-[baro_data[0,0]+0.8]*40
#t3=t[::96]

def fullboundary(vel,gps): #should return correct position in 1D. use boundary first to get vel from acc
    vel_right=[]#trimmed lists of vel to right of gps point for each gps point
    vel_left=[]#left of gps point
    n=int(vel.size/gps.size) #ratio between measurement freq. tells you how many vel to skip per gps
    
    for x in range (0,gps.size):
        vel_right.append([vel[(x*n):vel.size]])
        vel_left.append([vel[0:(x*n)]])
    right_fwd=[]#forward trajectory to right of gps with x0=x_gps for each gps point
    right_bkwd=[]#backward right
    left_fwd=[]
    left_bkwd=[]
    for x in range (0,gps.size):
        M=cumtrapz(vel_right[x],initial=gps[x])
        N=cumtrapz(vel_right[x][::-1],initial=0)
        O=cumtrapz(vel_left[x],initial=0)
        P=cumtrapz(vel_left[x][::-1],initial=gps[x])
        right_fwd.append(M)
        right_bkwd.append(N[::-1])
        left_fwd.append(O)
        left_bkwd.append(P[::-1])#gotta flip N and P cos the integrated data will be given bkwds to the fwd data 
    left_weight=[] # weighted trajectory of fwd and back for ONE gps point (to the right)
    left_weight_tot=[] #list of all weighted trajectories to right of gps points 
    right_weight=[]
    right_weight_tot=[]
     
    for t in range(0,gps.size):
        m=left_fwd[t].size #length of a specific left traj for a specific gps point
        k=right_fwd[t].size
        if m==1:
            left_weight.append(left_fwd[t][0][0])#do this to avoid dividing by zero in the weighting process; cos if m=1, m-1=0
        else:
            for x in range (0,m):
                left_f_w=((m-x-1)/(m-1))*left_fwd[t][0][x] #fwd weighted so most confident at points near gps point
                left_b_w=(x/(m-1))*left_bkwd[t][0][x] #bkwd weighted so most confident at end start of tot traj
                left_weight.append(np.mean([left_f_w,left_b_w])) #add each weighted point to the traj to make up one whole traj
        left_weight_tot.append(left_weight)#add specific weighted traj for specif gps point to set of all left traj
        if k==1:
            right_weight.append(right_fwd[t][0][0])
        else:
            for x in range (0,k):
                right_f_w=((k-x-1)/(k-1))*right_fwd[t][0][x]
                right_b_w=(x/(k-1))*right_bkwd[t][0][x]
                right_weight.append(np.mean([right_f_w,right_b_w]))  
                #print(np.mean([right_f_w,right_b_w]))
        right_weight_tot.append(right_weight)#same deal but for traj to the right of gps points 
        
    weight_tot=[]#list of the total weighted trajectories for each gps point
    for t in range (0,len(right_weight_tot)):
        weight_tot.append(left_weight_tot[t]+right_weight_tot[t])#just combine left and right trajectories
    full_trajectory=[]#weight each trajectory, with highest confidence for the points near the gps point
    for t in range(0,gps.size):
        for x in range(0,len(weight_tot)):
            if x<=round(n*t): # this is weighting for points to right of gps point
                w1=(len(weight_tot)-1-x)/(len(weight_tot)-round(n*t)-1) #stuff to the left of gps point
                full_trajectory.append(w1*weight_tot[t][x])                
            else:
                w2=x/(round(n*t)-1)
                full_trajectory.append(w2*weight_tot[t][x])
                
    return full_trajectory
#%%
VZ=boundary(AZ)
SZ=fullboundary((VZ+[66.19]*VZ.size),H)
plt.plot(SZ)