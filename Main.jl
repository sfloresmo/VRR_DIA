##################################################
include("Packages.jl") # Loading the packages
include("Args.jl") # Parsing the command-line
include("Constants.jl") # Physical constants
include("Elsasser.jl") # Elsasser coefficients
include("Coupling.jl") # Coupling coefficients
include("Prediction.jl") # Solving the prediction 
include("Order_1.jl") # Order 1 prediction corresponding to the DIA prediction
include("Order_2.jl") # Order 2 prediction corresponding to the next-order MSR diverging prediction