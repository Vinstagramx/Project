def veltopos(vel,pos,t): #should return correct position in 1D. use boundary first to get vel from acc
    vel_right=[]    #trimmed lists of vel to right of position point for each position point
    vel_left=[] #left of position point
    n = sp.round(len(vel)/len(pos)) #ratio between measurement freq. tells you how many vel to skip per position - has to be an integer
    
    for x in range (0,len(pos)):
        vel_right.extend([vel[(x*n):len(pos)]])
        vel_left.extend([vel[0:(x*n)]])
    right=[]#forward trajectory to right of position with x0=x_position for each a point
    right_flip=[]#backward right
    left=[]
    left_flip=[]
    for x in range (0,len(pos)):
        right.append(cumtrapz(vel_right[x],t,initial=pos[x]))

        right_flip_backwards = cumtrapz(vel_right[x][::-1],t,initial=0)
        right_flip.append(right_flip_backwards[::-1])

        left.append(cumtrapz(vel_left[x],t,initial=0))

        left_flip_backwards = cumtrapz(vel_left[x][::-1],t,initial=pos[x])
        left_flip.append(left_flip_backwards[::-1])     #gotta flip the backwards integrals cos the integrated data will be given bkwds to the fwd data 
   

    left_weight=[] # weighted trajectory of fwd and back for ONE altitude point (to the right)
    left_weight_tot=[] #list of all weighted trajectories to left of position points 
    right_weight=[]
    right_weight_tot=[]
     
    for t in range(0,len(pos)):
        m = len(left[t]) #length of a specific left traj for a specific position point
        k = len(right[t]).size

        if m==1:
            left_weight.append(left[t][0][0])#do this to avoid dividing by zero in the weighting process; cos if m=1, m-1=0
        else:
            for x in range (0,m):
                left_f_w=((m-x-1)/(m-1))*left[t][0][x] #fwd weighted so most confident at points near position point
                left_b_w=(x/(m-1))*left_flip[t][0][x] #bkwd weighted so most confident at end start of tot traj
                left_weight.append(np.mean([left_f_w,left_b_w])) #add each weighted point to the traj to make up one whole traj
        left_weight_tot.append(left_weight)#add specific weighted traj for specif position point to set of all left traj
        if k==1:
            right_weight.append(right[t][0][0])
        else:
            for x in range (0,k):
                right_f_w=((k-x-1)/(k-1))*right[t][0][x]
                right_b_w=(x/(k-1))*right_flip[t][0][x]
                right_weight.append(np.mean([right_f_w,right_b_w]))  
                #print(np.mean([right_f_w,right_b_w]))
        right_weight_tot.append(right_weight)#same deal but for traj to the right of position points 
        
    weight_tot=[]#list of the total weighted trajectories for each position point
    for t in range (0,len(right_weight_tot)):
        weight_tot.append(left_weight_tot[t]+right_weight_tot[t])#just combine left and right trajectories
    full_trajectory=[]#weight each trajectory, with highest confidence for the points near the position point
    for t in range(0,position.size):
        for x in range(0,len(weight_tot)):
            if x>=round(n*t): # this is weighting for points to right of position point
                w1=(len(weight_tot)-1-x)/(len(weight_tot)-round(n*t)-1) #stuff to the left of position point
                full_trajectory.append(w1*weight_tot[t][x])                
            else:
                w2=x/(round(n*t)-1)
                full_trajectory.append(w2*weight_tot[t][x])
                
    return full_trajectory
#%%