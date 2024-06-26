# julia --threads 8 Run_Prediction.jl --lmax 20 --coupling Quad --Npart 1000 --dt 1.0 --Nsteps 330 --order 1
##################################################
include("../src/Main.jl")
##################################################
println("LMAX | ",LMAX," | COUPLING | ",COUPLING," | ORDER | ",ORDER," | NPART | ",NPART)
println("DT | ",DT," | NSTEPS | ",NSTEPS," | T_C | ",T_C)
##################################################
# To dump the results
##################################################
function dump!()
    namefile = "data/data_Prediction_LMAX_"*string(LMAX)*"_COUPLING_"*COUPLING*"_ORDER_"*string(ORDER)*".hf5"
    #####
    file = h5open(namefile,"w") # Opening the file
    #####
    # Writing the parameters of the run
    write(file,"LMAX",LMAX)
    write(file,"COUPLING",COUPLING)
    write(file,"NPART",NPART)
    write(file,"DT",DT)
    write(file,"NSTEPS",NSTEPS)
    write(file,"ORDER",ORDER)
    write(file,"T_C",T_C)
    #####
    # Writing the data of the run
    write(file,"TAB_T",TAB_T)
    write(file,"TAB_R",TAB_R)
    #####
    close(file) # Output
end
##################################################
print("Computing |")
@time TAB_R!()
print("Dumping |")
@time dump!()
