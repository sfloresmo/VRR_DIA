using BenchmarkTools # To benchmark code precisely
using HDF5 # To use HDF5 files
using Polyester # Multi-threading with lower allocation imprint
using OffsetArrays # To have arrays that start at the index 0
using Memoize # To memoize the Elsasser coefficients
using WignerSymbols # To have access to Wigner symbols