using ArgParse
##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--lmax"
    help = "Maximum harmonic number"
    arg_type = Int64
    default = 10
    "--coupling"
    help = "Coupling considered: Quad/Asymp"
    arg_type = String
    default = "Quad"
    "--Npart"
    help = "Number of particles"
    arg_type = Int64
    default = 1000
    "--dt" #ATTENTION dt in physical units (NOT in T_C units)
    help = "Integration timestep"
    arg_type = Float64
    default = 1.0
    "--Nsteps"
    help = "Number of integration steps"
    arg_type = Int64
    default = 100
    "--order"
    help = "Order of the approximation: 1/2"
    arg_type = Int64
    default = 1
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################
const LMAX = parsed_args["lmax"] # Maximum harmonic number
const COUPLING = parsed_args["coupling"] # Coupling considered
const NPART = parsed_args["Npart"] # Number of particles
const DT = parsed_args["dt"] # Integration timestep #ATTENTION dt in physical units (NOT in T_C units)
const NSTEPS = parsed_args["Nsteps"] # Number of integration timesteps
const ORDER = parsed_args["order"] # Order of the approximation
##################################################
# Checking the sanity of some of the parameters
##################################################
if (COUPLING != "Quad") && (COUPLING != "Asymp")
    error("Supported COUPLING: Quad/Asymp")
end
#####
if (ORDER != 1) && (ORDER != 2)
    error("Supported ORDER: 1/2")
end