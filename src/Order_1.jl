##################################################
#Order 1 approximation corresponding to DIA prediction
#Computing the integral and the double sum appearing in Eq.(28)
##################################################
# Approximates the ORDER=1 convolution integral
# using the trapezoidal rule
# More precisely, we compute the integral
# INT[ R[l2,t0-t1] * R[l3,t0-t1] * R[l1,t1] ,{t1,0,t0}]
# This integral is used in the ORDER=1 prediction
##################################################
function calc_I_1(l1::Int64,l2::Int64,l3::Int64, n0::Int64)
    #####
    (n0 == 0) && return 0.0 # Returning 0.0 for an empty integral
    #####
    res = 0.0 # Initialising the result
    #####
    @inbounds for n1=0:n0 
        #####
        ((n1 == 0) || (n1 == n0)) ? eps_1 = 0.5 : eps_1 = 1.0 # Fixing the prefactor
        #####
        eps = eps_1
        #####
        res += eps * TAB_R[l2,n0-n1] * TAB_R[l3,n0-n1] * TAB_R[l1,n1]
        #####
    end
    #####
    res *= DT # Multiplying by the timestep
    #####
    return res # Output
end
##################################################
# Wrapped function that computes dR[l]/dt
# at the time t = n0 * DT
# Contribution from ORDER = 1
##################################################
function calc_dRdt_1(l1::Int64, n0::Int64)
    #####
    res = 0.0 # Initialising the result
    #####
    J1 = get_Jl(l1)
    #####
    for l2=0:LMAX
        #####
        J2 = get_Jl(l2)
        (J1 == J2) && continue
        #####
        for l3=0:LMAX
            #####
            J3 = get_Jl(l3)
            (J2 == J3) && continue
            #####
            EL_123 = get_ElsasserL(l1,l2,l3)
            (EL_123 == 0.0) && continue
            #####
            # Computing the prefactor in front of the integral,
            # but within the sum over harmonics
            pref = (J2 - J3) * (J2 - J1) * EL_123 * EL_123
            #####
            # Adding the contribution from the coupling (l1,l2,l3)
            # In practice, we already tested for all the possible exclusions,
            # i.e., normally, we are sure that pref is non-zero,
            # i.e., we are sure that it is worth computing this integral
            res += pref * calc_I_1(l1,l2,l3,n0)
        end
    end
    #####
    # Multiplying by the appropriate prefactor
    pref_out = - (1.0/(2.0*l1+1.0))*(NPART/(4.0*pi)) # ATTENTION, not to forget the minus sign
    res *= pref_out
    #####
    return res # Output
end