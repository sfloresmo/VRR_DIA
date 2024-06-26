Vector Resonant Relaxation (VRR)
-
Stars orbiting a supermassive black hole
in the center of galaxies undergo very efficient
diffusion in their orbital orientations. This is "Vector Resonant Relaxation". 
We present a code that computes an approximate prediction of the two-point correlation function. The first-order approximation (order 1) corresponds to the so-called Direct Interaction Approxiamtion (DIA). 
The next-order (order 2) diverging prediction is also computed in this code.

Packages
-
This code is written in Julia. To install Julia: https://julialang.org/downloads/platform/

List of Julia packages required to run the code:
- BenchmarkTools 
- HDF5 
- Polyester 
- OffsetArrays 
- Memoize
- WignerSymbols 
- Plots

Run
-
First, open the run folder in the terminal:
```sh
cd run
```
Second, generate the data of the prediction by writing the command (computing time < 2s): 

```sh
julia --threads 8 Run_Prediction.jl --lmax 20 --coupling Quad --Npart 1000 --dt 1.0 --Nsteps 330 --order 1
```
Third, generate the figure from the data write:
```sh
julia Run_Create_Figs.jl --lmax 20 --coupling Quad --Npart 1000 --dt 1.0 --Nsteps 330 --order 1
```
You can generate the data/figure for couplings Quad/Asymp and orders 1/2 of approximation.

Figure
-
![figure_Prediction_LMAX_20_COUPLING_Asymp_ORDER_1](https://github.com/sfloresmo/VRR_DIA/assets/172500564/32f48f07-8e31-44aa-bc8e-8bf1c4ba9510)
