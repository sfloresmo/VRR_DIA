##################################################
# In this code, we solve the prediction equations (28) and (G2)
# within the Fluctuation-Dissipation Theorem 
##################################################
# Arrays used for the integration
# ATTENTION, the arrays are off-set
const TAB_T    = OffsetArray(zeros(Float64,NSTEPS+1),0:NSTEPS) # Table of the considered time. 
const TAB_R    = OffsetArray(zeros(Float64,LMAX+1,NSTEPS+1),0:LMAX,0:NSTEPS) # Values of the response functions
#####
const TAB_DRDT = OffsetArray(zeros(Float64,LMAX+1),0:LMAX) # Container to keep the value of the gradients dR/dt
##################################################
# Initializing the times
##################################################
function TAB_T!()
    for n=0:NSTEPS
        TAB_T[n] = n * DT
    end
end
##################################################
TAB_T!() # Filling in the values of the times
##################################################
# Wrapped function that computes dR/dt
# at the appropriate order
##################################################
function calc_dRdt_wrapped_1(l::Int64,n0::Int64)
    return calc_dRdt_1(l,n0)
end
##################################################
function calc_dRdt_wrapped_2(l::Int64,n0::Int64)
    return calc_dRdt_1(l,n0) + calc_dRdt_2(l,n0)
end
##################################################
# Picking the appropriate function
# for the considered ORDER
# ATTENTION, not to forget the `const'
##################################################
if     (ORDER == 1)
    const calc_dRdt = calc_dRdt_wrapped_1
elseif (ORDER == 2)
    const calc_dRdt = calc_dRdt_wrapped_2
end
##################################################
# Updates for one timestep assuming that we know
# the value of R(l,n*DT) for 0<=n<=n0
# Phrased differently, it computes R[l,(n0+1)*DT]
##################################################
function update_TAB_R!(n0::Int64)
    #####
    @batch for l=0:LMAX
        #####
        dRdt_n = calc_dRdt(l,n0) # Computing the value of dR[l,n0]/dt
        TAB_DRDT[l] = dRdt_n # Keeping the value as it will be used later on
        #####
        R_n = TAB_R[l,n0] # Value of R[l,n0]
        hR_np1 = R_n + DT * dRdt_n # Predicting the value of R[l,n0+1]
        #####
        TAB_R[l,n0+1] = hR_np1 # Filling in with the predicted value
    end
    #####
    @batch for l=0:LMAX
        #####
        dRdt_np1 = calc_dRdt(l,n0+1) # Computing the value of dR[l,n+1]/dt
        #####
        R_n = TAB_R[l,n0] # Value of R[l,n]
        dRdt_n = TAB_DRDT[l] # Stored value of dR[l,n]/dt
        #####
        R_np1 = R_n + 0.5 * DT * (dRdt_n + dRdt_np1) # Correcting the value of R[l,n+1]
        #####
        TAB_R[l,n0+1] = R_np1 # Updating the array
    end
end
##################################################
# Initializing the values of the response functions
##################################################
function init_TAB_R!()
    for l=0:LMAX
        TAB_R[l,0] = 1.0 # The response functions are all equal to 1 for t=0
    end
end
##################################################
# Computing all the values of TAB_R
##################################################
function TAB_R!()
    #####
    init_TAB_R!() # Initialising TAB_R
    #####
    for n0=0:(NSTEPS-1) # Performing the loop
        update_TAB_R!(n0) # Computing TAB_R[:,n0+1]
         println("n = ",n0,"/",NSTEPS," | t/Tc = ", n0*DT/T_C)
         println("**************")
    end
    #####
end