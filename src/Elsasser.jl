##################################################
# Implementing the isotropic part of the Elsasser coefficients
# We follow Appendix in Fouvry+2019
##################################################
# Coefficient Lambda appearing in E^L
# See Eq. (56) in Fouvry+2019
##################################################
function get_Lambda_Elsasser(la::Int64,lb::Int64,lc::Int64)
    return sqrt(((2.0*la+1.0)*(2.0*lb+1.0)*(2.0*lc+1.0))/(4.0 * pi))
end
##################################################
# Coefficient Delta appearing in E^L
# See Eq. (56) in Fouvry+2019
##################################################
function get_Delta_Elsasser(la::Int64,lb::Int64,lc::Int64)
    return sqrt(((la+lb+lc+2.0)*(la+lb+lc+4.0))/(4.0*(la+lb+lc+3.0))) *
           sqrt((la+lb-lc+1.0)*(lc+la-lb+1.0)*(lb+lc-la+1.0))
end
##################################################
# Wigner 3j symbol appearing in E^L
# See Eq. (55) in Fouvry+2019
##################################################
function get_Wigner3j_Elsasser(la::Int64,lb::Int64,lc::Int64)
    return wigner3j(la+1,lb+1,lc+1,0,0,0)
end
##################################################
# Returns the E^L coefficients, E^L[la,lb,lc]
# ATTENTION, the function is memoised
# to be rapidly accessed afterwards
##################################################
@memoize function get_ElsasserL(la::Int64,lb::Int64,lc::Int64)
    #####
    # Applying some simple exclusion rules
    # The sum of harmonics has to be odd,
    # for the coefficients to have a chance of being non-zero
    if (iseven(la+lb+lc))
        return 0.0
    end
    #####
    # The coefficients have to comply with a strict triangular inequality
    # for the coefficients to have a chance to be non-zero
    if (!(abs(la-lb) < lc < (la+lb)))
        return 0.0
    end
    #####
    # Following Eq. (55) of Fouvry+2019,
    # the Elsasser coefficients read
    # E^L = Lambda * Delta * Wigner
    val_Lambda = get_Lambda_Elsasser(la,lb,lc)
    val_Delta  = get_Delta_Elsasser(la,lb,lc)
    val_Wigner = get_Wigner3j_Elsasser(la,lb,lc)
    #####
    val_EL = val_Lambda * val_Delta * val_Wigner # Elsasser coefficients
    #####
    return val_EL # Output
end
##################################################
# Returns the W coefficients, W[la,lc,ld,lg,le,lf], given by
# = (-1)^(lc+le) * wigner6j[la,lc,ld,lg,le,lf]
# ATTENTION, the function is memoised to be rapidly accessed afterwards
##################################################
@memoize function get_W_Elsasser(la::Int64,lc::Int64,ld::Int64,
                                 lg::Int64,le::Int64,lf::Int64)
    #####
    res = (-1)^(lc+le) *  wigner6j(la,lc,ld,lg,le,lf) # Computing the exact value
    res = convert(Float64,res) # Converting to Float64
    #####
    return res # Output
end