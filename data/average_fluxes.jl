using JLD2
using Oceananigans
using OceanScalingTests
using OceanScalingTests: exponential_z_faces
using SeawaterPolynomials
using SeawaterPolynomials: TEOS10EquationOfState

using Oceananigans.Fields: ConstantField
using Oceananigans.Models: seawater_density
using Oceananigans.AbstractOperations: GridMetricOperation
using CUDA 

Na = (4320, 1800, 1)
Ny = (4320, 1801, 1)

grid = LatitudeLongitudeGrid(CPU(), Float32; size = Na, latitude = (-75, 75), longitude = (-180, 180), z = (0, 1))

u = XFaceField(grid)
v = YFaceField(grid) 
T = CenterField(grid) 
S = CenterField(grid) 

Q = CenterField(grid)
F = CenterField(grid)
τu = XFaceField(grid)
τv = YFaceField(grid)

myfiles = readdir("./")
myfiles = filter(x -> length(x) >= 7, myfiles)
myfiles = filter(x -> x[1:7] == "fluxes_", myfiles)

@inline remove_last_character(s) = s[1:end-1]
myfiles = remove_last_character.(myfiles)

numbers = parse.(Int, filter.(isdigit, myfiles))
numbers = sort(numbers)

for num in numbers[1:end]
   GC.gc()
   @info "doing file number $num"
   file = jldopen("fluxes_" * string(num) * ".jld2")
   for n in 1:5
      set!(u, file["τx"][:, :, n])
      set!(v, file["τy"][:, :, n])
      set!(T, file["Qs"][:, :, n])
      set!(S, file["Fs"][:, :, n])

      τu .+= u / 5 / 79
      τv .+= v / 5 / 79
      Q  .+= T / 5 / 79
      F  .+= S / 5 / 79
   end
end

jldsave("average_fluxes.jld2", τu=τu, τv=τv, Q=Q, F=F)
