##################################################
# Coupling coefficient in the case of the quadrupolar model
# In that case, only the \ell=2 coupling term is non-zero
##################################################
function get_Jl_quad(l::Int64)
    (l == 2) ? J2_COUPLING : return 0.0
end
##################################################
# Coupling coefficient for arbitrary \ell
# Following Kocsis&Tremaine(2015) or Takacs&Kocsis(2017),
# with the asymptotic scaling
# J_\ell = J_2 / (\ell/2)^2
# with J_2 given by the \ell=2 coupling coefficient
# in the case of a single population system
##################################################
function get_Jl_asymp(l::Int64)
    isodd(l) && return 0.0 # The interaction is zero if \ell is off
    (l == 0) && return 0.0 # The \ell=0 interaction coefficient can be taken to 0
    return J2_COUPLING/((l/2)^2) # Returning the coefficient with its asymptotic scaling
end
##################################################
# Picking the appropriate coupling coefficient
# depending on the interaction model that is considered
##################################################
if     (COUPLING == "Quad")
    const get_Jl = get_Jl_quad
else (COUPLING == "Asymp")
    const get_Jl = get_Jl_asymp
end