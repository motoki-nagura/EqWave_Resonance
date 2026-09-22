#=======================================================================
      Prepare forcing data for idealized numerical model experiments
========================================================================#

#================================================================================
                                  Parameters
================================================================================#

#const year        = 1       #   Year.  Year 1 == 2000,  Year 2 = 2001, etc.
const year = parse(Int, ARGS[1])

#const wind_period = 180     #   Period [days]. 180, 90 or 135
#const wind_period =  90
#const wind_period =  135
const wind_period = parse(Int, ARGS[2])

println(" =====>  ARGS input from loop.sh")
println("year = ",year,",  wind_period = ",wind_period)
println()


#     For Wind forcing


#     Grids
const imt, jmt, km = 622, 402, 100
file_name_grid     = "/S/data01/G4006/y0330/data/ofes_exp/settings/box_bounded/topog/grid.dta.out"
file_name_grid_lib = "/S/data01/G4006/y0330/code/misc/ofes_exp/ReadGrids/ForJulia/libmysub_box_bounded.so"

#     Read in timesteps
infile_time = "/S/data01/G4006/y0330/data/misc/jra55-do/dayavg/uas/"*
   "uas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-4-0_gr_201001010000-201012312100.day.nc"

#================================================================================
                           Read grids from data file   {{{
================================================================================#
println(" ====>  Read in grids from data file")
println("IN <<< "*file_name_grid)

isfile("./grid.dta.out") && rm("./grid.dta.out")
symlink(file_name_grid, "./grid.dta.out")

xt, yt = zeros(Float64, imt), zeros(Float64, jmt)
xu, yu = zeros(Float64, imt), zeros(Float64, jmt)
zt, zw = zeros(Float64, km),  zeros(Float64, km)

ccall((:read_grid_dta_, file_name_grid_lib),
       Cvoid,
       (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
       xt,yt,xu,yu,zt,zw)

lon, lat = xu, yu
println("lon = ",lon[1],", ",lon[2],", ..., ",lon[end])

println() #}}}
#================================================================================
                           Read in timesteps  {{{
================================================================================#
println(" ====>  Read in timesteps")
using NCDatasets

println("IN <<< "*infile_time)

ds = NCDataset(infile_time,"r")
time = Array(ds["time"])
close(ds)

time = time .- time[1]                    #  if input data is daily
#time = (time .- time[1]) ./ 24.          #                   hourly

time = getfield.(time, :value)
time = time / 86_400_000                  #  Milli seconds --> days

time = time .+ (year-1) *365.

println("  time = ",time[1],", ",time[2],", ... , ",time[end],", ",length(time)," data points")

println() #}}}
#================================================================================
                            Set τˣ  {{{
================================================================================#
println(" ====>  Set τˣ")

if      wind_period == 180
    global Amp, lon_cent, lon_width, lat_width = 0.15, 30., 50., 15.
elseif  wind_period == 90
    global Amp, lon_cent, lon_width, lat_width = 0.06, 45., 40., 15.
else
    global Amp, lon_cent, lon_width, lat_width = 0.15, 30., 50., 15.
   #   For other periods, the spatial structure of 180-dy winds is used.
end

#------   Horizontal structure

x = 2π .* ( xu .- lon_cent ) ./lon_width
xfunc_model = 0.5 .*( cos.(x) .+ 1. )
xfunc_model[ x .< -π ] .= 0.
xfunc_model[ x .>  π ] .= 0.

yfunc_model = exp.( -0.5 .* (yu .^2) ./lat_width )

func_model = xfunc_model * yfunc_model'     #  matrix multiplication
func_model = Amp .* func_model

#------   Temporal function

t_prof = sin.( 2π .*time ./wind_period )

τˣ = func_model[:, :, :] .* reshape(t_prof, 1, 1, :)
τʸ = fill(0., size(τˣ))

println()  #  }}}
#================================================================================
                            Write out to NetCDF files {{{
================================================================================#
println(" ====>  Write out to NetCDF files")
using NetCDF
using Printf

year_actual = year + 2000 - 1

#-----

ext = "_" *string(wind_period) *"dy" *"_" *string(year_actual, base = 10, pad = 4)

file_name = "./netcdfs/taux" * ext * ".nc"
println("  OUT >>> " * file_name)
isfile(file_name) && rm(file_name)
ds = NCDataset(file_name,"c")
v = defVar(ds,"taux", τˣ,   ("lon","lat","time"))
x = defVar(ds,"lon",  xu,   ("lon",))
y = defVar(ds,"lat",  yu,   ("lat",))
t = defVar(ds,"time", time, ("time",))
v.attrib["longname"]   = "Zonal wind stress"
v.attrib["units"   ]   = "dyn cm^-2"
x.attrib["longname"]   = "Longitude"
x.attrib["units"   ]   = "degrees east"
y.attrib["longname"]   = "Latitude"
y.attrib["units"   ]   = "degrees north"
t.attrib["longname"]   = "Time"
t.attrib["units"   ]   = "days since 2000-1-1 00:00:00"
close(ds)

file_name = "./netcdfs/tauy" * ext * ".nc"
println("  OUT >>> " * file_name)
isfile(file_name) && rm(file_name)
ds = NCDataset(file_name,"c")
v = defVar(ds,"tauy", τʸ,   ("lon","lat","time"))
x = defVar(ds,"lon",  xu,   ("lon",))
y = defVar(ds,"lat",  yu,   ("lat",))
t = defVar(ds,"time", time, ("time",))
v.attrib["longname"]   = "Meridional wind stress"
v.attrib["units"   ]   = "dyn cm^-2"
x.attrib["longname"]   = "Longitude"
x.attrib["units"   ]   = "degrees east"
y.attrib["longname"]   = "Latitude"
y.attrib["units"   ]   = "degrees north"
t.attrib["longname"]   = "Time"
t.attrib["units"   ]   = "days since 0001-1-1 00:00:00"
close(ds)

println() #  }}}

