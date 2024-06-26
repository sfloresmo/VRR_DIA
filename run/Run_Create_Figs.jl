# julia Run_Create_Figs.jl --lmax 20 --coupling Quad --Npart 1000 --dt 1.0 --Nsteps 330 --order 1
##################################################
include("../src/Main.jl")
##################################################
println("LMAX | ",LMAX," | COUPLING | ",COUPLING," | ORDER | ",ORDER," | NPART | ",NPART)
println("DT | ",DT," | NSTEPS | ",NSTEPS," | T_C | ",T_C)
##################################################
# Creating the plot
plot(time2, values2, label="l=2", xlabel="t/Tc", ylabel="R_l")
plot!(time3, values3, label="l=3")
plot!(time4, values4, label="l=4")

# Save the plot 
figs_directory = string("../figs/")
savefig(figs_directory * "figure_Prediction_LMAX_"*string(LMAX)*"_COUPLING_"*COUPLING*"_ORDER_"*string(ORDER)*".png")