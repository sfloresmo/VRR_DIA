function read_prediction(lmax::Int64, coupling::String, order::Int64, l::Int64)
    # Name of the file
    namefile = "../data/data_Prediction_LMAX_"*string(lmax)*"_COUPLING_"*coupling*"_ORDER_"*string(order)*".hf5"
    
    # Coherence time
    Tc = h5read(namefile, "T_C")
    
    # Time
    tabt = h5read(namefile, "TAB_T")
    
    # Rescaling the time with respect to Tc
    tabt = tabt / Tc
    
    # Reading the response functions
    # It has the shape {t, 0 <= l <= LMAX}
    tabR = h5read(namefile, "TAB_R")
    
    # Extracting the appropriate harmonics
    # ATTENTION, the harmonics start at l = 0
    tabRl = tabR[l+1, :]
    
    # Output
    return tabt, tabRl
end

# Reading the prediction for the harmonics l=2, 3, 4
time2, values2 = read_prediction(LMAX, COUPLING, ORDER, 2)
time3, values3 = read_prediction(LMAX, COUPLING, ORDER, 3)
time4, values4 = read_prediction(LMAX, COUPLING, ORDER, 4)




