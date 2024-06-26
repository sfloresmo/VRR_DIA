##################################################
# Order 2 approximation corresponding to next-order MSR diverging prediction
# Computing the triple integral and the five sums appearing in the second term of Eq.(G2)
#############################################
# We compute the integral
# INT[ R[l1,t3] * R[l2,t0-t2] * R[l3,t0-t1] * R[l4,t1-t2] * R[l5,t1-t3] * R[l6,t2-t3] 
# * H(-) * H(-) * H(-) , {t1,0,t0},{t2,0,t0},{t3,0,t0}]
# for each of the four terms in the triple integral,
# with H the usual Heaviside function
##########################################################
function get_R(l::Int64, n::Int64)# abs() prevents problem with index in TAB_R
    return TAB_R[l, abs(n)]
end
##########################
# Computing each one of the four terms appearing in the triple integral 
# without the prefactors
# Terms 1,2,3,4 numbered in the order that they appear in Eq.(G2)
##########################
# Term 1
##########################
function term_1(l1::Int64,l2::Int64,l3::Int64,l4::Int64,l5::Int64,l6::Int64,n0::Int64)
    #####
    (n0 == 0) && return 0.0 # Returning 0.0 for an empty integral
    #####
    res = 0.0 # Initializing the result
    #####
    for n1=0:n0 
    #####
        ((n1 == 0) || (n1 == n0)) ? eps_1 = 0.5 : eps_1 = 1.0 
        #####
        R_l3 = get_R(l3,n0-n1)
            ####
            for n2=1:n0 
                #####
                (n2 == n0) ? eps_2 = 0.5 : eps_2 = 1.0 
                #####
                R_l2 = get_R(l2,n0-n2)
                R_l4 = get_R(l4,n1-n2) 
                ####
                for n3=0:n2 
                    #####
                    ((n3 == 0) || (n3 == n2)) ? eps_3 = 0.5 : eps_3 = 1.0 
                    #####
                    R_l1 = get_R(l1,n3) 
                    R_l5 = get_R(l5,n1-n3) 
                    R_l6 = get_R(l6,n2-n3) 
                    res += eps_1 * eps_2 * eps_3 * R_l1 * R_l2 * R_l3 * R_l4 * R_l5 * R_l6
                    #####
                end
            end
    end
    #####
    res *= DT * DT * DT # Multiplying by the timestep cubed
    #####
    return res # Output
end
#############################
# Term 2
#############################
function term_2(l1::Int64,l2::Int64,l3::Int64,l4::Int64,l5::Int64,l6::Int64,n0::Int64)
        #####
        (n0 == 0) && return 0.0 # Returning 0.0 for an empty integral
        #####
        res = 0.0 # Initializing the result
        #####
        for n2=1:n0 
        #####
            (n2 == n0) ? eps_2 = 0.5 : eps_2 = 1.0 
            #####
            R_l2 = get_R(l2,n0-n2)
                ####
                for n3=1:n2
                    #####
                    (n3 == n2) ? eps_3 = 0.5 : eps_3 = 1.0 
                    #####
                    R_l1 = get_R(l1,n3) 
                    R_l6 = get_R(l6,n2-n3)
                    ####
                    for n1=0:n3 
                        #####
                        ((n1 == 0) || (n1 == n3)) ? eps_1 = 0.5 : eps_1 = 1.0 
                        #####
                        R_l3 = get_R(l3,n0-n1)
                        R_l4 = get_R(l4,n1-n2)
                        R_l5 = get_R(l5,n1-n3) 
                        res += eps_1 * eps_2 * eps_3 * R_l1 * R_l2 * R_l3 * R_l4 * R_l5 * R_l6
                        #####
                    end
                end
        end
        #####
        res *= DT * DT * DT # Multiplying by the timestep cubed
        #####
    #####
    return res # Output
end
#############################
# Term 3
#############################
function term_3(l1::Int64,l2::Int64,l3::Int64,l4::Int64,l5::Int64,l6::Int64,n0::Int64)
    #####
    (n0 == 0) && return 0.0 # Returning 0.0 for an empty integral
    #####
    res = 0.0 # Initializing the result
    #####
    for n1=1:n0 
    #####
        (n1 == n0) ? eps_1 = 0.5 : eps_1 = 1.0 
        #####
        R_l3 = get_R(l3,n0-n1)
            ####
            for n2=1:n1
                #####
                (n2 == n1) ? eps_2 = 0.5 : eps_2 = 1.0 
                #####
                R_l2 = get_R(l2,n0-n2)
                R_l4 = get_R(l4,n1-n2)
                ####
                for n3=0:n2 
                    #####
                    ((n3 == 0) || (n3 == n2)) ? eps_3 = 0.5 : eps_3 = 1.0 
                    #####
                    R_l1 = get_R(l1,n3) 
                    R_l5 = get_R(l5,n1-n3) 
                    R_l6 = get_R(l6,n2-n3) 
                    res += eps_1 * eps_2 * eps_3 * R_l1 * R_l2 * R_l3 * R_l4 * R_l5 * R_l6
                    #####
                end
            end
    end
    #####
    res *= DT * DT * DT # Multiplying by the timestep cubed
    #####
    return res # Output
end
#############################
# Term 4
#############################
function term_4(l1::Int64,l2::Int64,l3::Int64,l4::Int64,l5::Int64,l6::Int64,n0::Int64)
   #####
   (n0 == 0) && return 0.0 # Returning 0.0 for an empty integral
   #####
   res = 0.0 # Initializing the result
   #####
    for n2=1:n0 
    #####
       (n2 == n0) ? eps_2 = 0.5 : eps_2 = 1.0 
       #####
       R_l2 = get_R(l2,n0-n2)
           ####
           for n3=0:n2 
               #####
               ((n3 == 0) ||(n3 == n2)) ? eps_3 = 0.5 : eps_3 = 1.0 
               #####
               R_l1 = get_R(l1,n3) 
               R_l6 = get_R(l6,n2-n3) 
               ####
               for n1=0:n2 
                   #####
                   ((n1 == 0) || (n1 == n2)) ? eps_1 = 0.5 : eps_1 = 1.0 
                   #####
                   R_l5 = get_R(l5,n1-n3) 
                   R_l4 = get_R(l4,n1-n2) 
                   R_l3 = get_R(l3,n0-n1)
                   res += eps_1 * eps_2 * eps_3 * R_l1 * R_l2 * R_l3 * R_l4 * R_l5 * R_l6
                   #####
               end
           end
    end
    #####
    res *= DT * DT * DT # Multiplying by the timestep cubed
    #####
    return res # Output
end
####################
#Summing the four terms with their corresponding prefactors
####################
function calc_I_2(l1::Int64,l2::Int64,l3::Int64,l4::Int64,l5::Int64,l6::Int64,n0::Int64)
    J2 = get_Jl(l2)
    J3 = get_Jl(l3)
    J4 = get_Jl(l4)
    J5 = get_Jl(l5)
    J6 = get_Jl(l6)

    pref_1 = (J5 - J4) * (J4 - J6) # prefactor of term 1
    pref_2 = (J4 - J3) * (J4 - J6) # prefactor of term 2
    pref_3 = (J5 - J4) * (J6 - J2) # prefactor of term 3
    pref_4 = (J4 - J6) * (J3 - J5) # prefactor of term 4

    (pref_1 == 0.0) ? int_1 = 0.0 : int_1 = pref_1 * term_1(l1,l2,l3,l4,l5,l6,n0)
    (pref_2 == 0.0) ? int_2 = 0.0 : int_2 = pref_2 * term_2(l1,l2,l3,l4,l5,l6,n0)
    (pref_3 == 0.0) ? int_3 = 0.0 : int_3 = pref_3 * term_3(l1,l2,l3,l4,l5,l6,n0)
    (pref_4 == 0.0) ? int_4 = 0.0 : int_4 = pref_4 * term_4(l1,l2,l3,l4,l5,l6,n0)

    return int_1 + int_2 + int_3 + int_4
end
##################################################
# Wrapped function that computes dR[l1]/dt0
# at the time t0 = n0 * DT
# Contribution from ORDER = 2, second term of Eq. (G2)
##################################################
function calc_dRdt_2(l1::Int64, n0::Int64)
    #####
    res = 0.0 # Initializing the result
    #####
    J1 = get_Jl(l1)
    #####
    for l5=0:LMAX
        #####
        J5 = get_Jl(l5)
        (J5 == J1) && continue
        #####
        for l2=0:LMAX
            #####
            J2 = get_Jl(l2)
            #####
            for l3=0:LMAX
                J3 = get_Jl(l3)
                (J2 == J3) && continue
                #####
                EL_123 = get_ElsasserL(l1,l2,l3)
                (EL_123 == 0.0) && continue
                #####
                for l4=0:LMAX
                    #####
                    EL_453 = get_ElsasserL(l4,l5,l3)
                    (EL_453 == 0.0) && continue
                    #####
                    for l6=0:LMAX
                        #####
                        EL_156 = get_ElsasserL(l1,l5,l6)
                        (EL_156 == 0.0) && continue
                        #####
                        EL_426 = get_ElsasserL(l4,l2,l6)
                        (EL_426 == 0.0) && continue
                        #####
                        W_123456 = get_W_Elsasser(l1,l2,l3,l4,l5,l6) # Result of the contraction
                        (W_123456 == 0.0) && continue
                        #####
                        # Computing the prefactor in front of the triple integral,
                        # but within the sum over harmonics
                        #####
                        pref = (J1 - J5) * (J2 - J3) * EL_123 * EL_156 * EL_426 * EL_453 * W_123456
                        #####
                        # In practice, we already tested for all the possible exclusions,
                        # i.e., normally, we are sure that pref is non-zero,
                        # i.e., we are sure that it is worth computing this integral
                        res += pref * calc_I_2(l1,l2,l3,l4,l5,l6,n0)
                    end
                end
            end
        end
    end
    #####
    # Multiplying by the appropriate outside prefactor
    pref_out = (1.0/(2.0*l1+1.0)) * (NPART/(4.0*pi))^(2)
    res *= pref_out
    #####
    return res # Output
end