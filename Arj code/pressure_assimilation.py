import numpy as np

AZ=fftfilter(t,accelreal[:,0][5240:9079])[3]#trimmed data we want to work with
H=np.asarray(h[55:95])
#
#t2=baro_data[:,0][55:95]-[baro_data[0,0]+0.8]*40
#t3=t[::96]

def fullboundary(vel,gps): #should return correct position in 1D. use boundary first to get vel from acc
    vel_right=[]#trimmed lists of vel to right of gps point for each gps point
    vel_left=[]#left of gps point
    n=vel.size/gps.size #ratio between measurement freq. tells you how many vel to skip per gps
    
    for x in range (0,gps.size):
        vel_right.append([vel[round(x*n):vel.size]])
        vel_left.append([vel[0:round(x*n)]])
    right_fwd=[]#forward trajectory to right of gps with x0=x_gps for each gps point
    right_bkwd=[]#backward right
    left_fwd=[]
    left_bkwd=[]
    for x in range (0,gps.size):
        M=cumtrapz(vel_right[x],initial=gps[x])
        N=cumtrapz(vel_right[::-1][x],initial=0)
        O=cumtrapz(vel_left[x],initial=0)
        P=cumtrapz(vel_left[::-1][x],initial=gps[x])
        right_fwd.append(M)
        right_bkwd.append(N[::-1])
        left_fwd.append(O)
        left_bkwd.append(P[::-1])
    left_weight=[] # weighted trajectory of fwd and back for ONE gps point (to the right)
    left_weight_tot=[] #list of all weighted trajectories to right of gps points 
    right_weight=[]
    right_weight_tot=[]
     
    for t in range(0,gps.size):
        m=left_fwd[t,:].size #length of a specific left traj for a specific gps point
        k=right_fwd[t,:].size
        for x in range (0,m):
            left_f_w=((m-x-1)/(m-1))*left_fwd[t,x] #fwd weighted so most confident at points near gps point
            left_b_w=(x/(m-1))*left_bkwd[t,x] #bkwd weighted so most confident at end start of tot traj
            left_weight.append(np.mean([left_f_w,left_b_w])) #add each weighted point to the traj to make up one whole traj
        left_weight_tot.append(left_weight)#add specific weighted traj for specif gps point to set of all left traj
        for x in range (0,k):
            right_f_w=((k-x-1)/(k-1))*right_fwd[t,x]
            right_b_w=(x/(k-1))*right_bkwd[t,x]
            right_weight.append(np.mean([right_f_w,right_b_w]))  
        right_weight_tot.append(right_weight)#same deal but for traj to the right of gps points 
        
    weight_tot=[]#list of the total weighted trajectories for each gps point
    for t in range (0,right_weight_tot):
        weight_tot.append(left_weight_tot[t]+right_weight_tot[t])
    full_trajectory=[]
    for t in range(0,gps.size):
        for x in range(0,weight_tot.size):
            if x>round(n*t):
                w1=x/(round(n*t)-1)
                full_trajectory.append(w1*weight_tot[x])
            else:
                w2=(weight_tot.size-1-x)/(weight_tot.size-round(n*t)-1)
                full_trajectory.append(w2*weight_tot[x])
    return full_trajectory
    #%%            
VZ=boundary(AZ)
#SZ=fullboundary(VZ,h)

#plt.plot(SZ)
