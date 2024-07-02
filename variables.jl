using Plots
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames
using Dates
using AstroLib
using SatelliteToolbox                                     # LIV 

#---Constants---#
const Isp = 800 # [s] Specific impulse
const Tmax = 0.0022 # [N] Maximun thrust
const μ = 3.986e5 #Orbital parameter [km^3 / s^2]
const g0 = 9.81 # [m/s^2] Standar gravity
const Re = 6378 # [Km] Earth Radio
const J2 = 1082.639e-6
const sigma = 5.67e-8 # Stefan Boltzmann's constant
const C_d = 2                  # Drag coefficient -->            LIV
const Q_sun = 1354 # [W/m^2] Sun Heatflux
const af = 0.19 #Albedo constant
const Q_ir = 228 # Earth Heatflux [W/m^2]

#---Propagator setup---#

Δt = 60 #Step size desired
Days = 20 #Simulation days
Tsim = 86400*Days # Total simulation time in seconds
N = round(Int,(Tsim/(Δt))+1) #number of points required
J2_per = true #J2 perturbation switch on/off
EclipseC = false # Thrust coasting due eclipse
Drag_per = true              # Drag perturbation switch on/off --->          LIV


#---Initial conditions---#
initial = [7071, 0.02, deg2rad(97.7), deg2rad(50.0), deg2rad(50.0), deg2rad(0.0), 24] # Initial orbit [a,e,i,ω,Ω,θ,m]
target = [7081, 0.01, deg2rad(97.75), 0.0, 0.0, 0.0] # target elements [a,e,i]
weights = [4.0, 1.0, 1.0, 0.0, 0.0] #Qlaw weights vector
Temp_init = 313.15 #Kelvin

#---Mission Start Date--#
Initial_Day = 20
Initial_Month = 4
Initial_Year = 2000

#---Arrays---#
n = 8 # Number of nodes for thermal
t = LinRange(0, Tsim, N) # time vector
x_mee = zeros((N, length(initial))) # mee array
x_coe = zeros((N, length(initial))) # coe array
controls = zeros((N,3)) # controls array
Temp=zeros(length(t),n) #temperature array
heat=zeros(length(t),n) #heat array
Earth_2_normal = zeros(n) 
Sun_2_normal = zeros(n)
Sun_2_Earth = zeros(n)


#------Material Thermal properties-----#
alpha =   [0.96, 0.96, 0.96, 0.03, 0.96, 0.96, 0.03] #absorptivity
epsilon = [0.85, 0.85, 0.85, 0.95, 0.25, 0.25, 0.95] #emissivity
cp = 961 #[J/kg k] Specific Heat
k = 400 #Contact conductance [W/k]
rho = 2770 #Material density [Kg/m^3]

#---Thruster thermal properties---#
kₜ = 400 #[W/m K] Copper thermal conductance
Aₜ = 0.1*0.1 #[m^2] Thruster conductance area
Δxₜ = 0.01 # [m] Thruster's mounting plate thickness

#-----Node's geometry------#
D=0.2  #depth [m]
L=0.35 #length [m]
W=0.01 # wall's Width [m]
Areas=[D*L, D*L, D*D, D*D, D*L, D*L, Aₜ, Aₜ] #Area's array [m^2]
C = rho.*Areas.*W.*cp  #Constant [J/K]
C[8] = 14*900

#---Thermal restrictions---#
Temp_sup = 423.15 # [K] Max temperature for thruster
Temp_inf = 263.15 # [K] Min temperature for thruster

Temp_max = 323.15 # [K] Max temperature for Res+PPU
Temp_min = 283.15 # [K] Min temperature for Res+PPU

#---Thruster power---#
η = 0.173 # Total efficiency
ηₚ = 0.70 # PPU efficiency
P_in = 50 # [W] Input power
P_out_ppu = ηₚ*P_in # [W] PPU power outpu
Pd = η*P_in # [W] Discharge power

#-------------- LIV:

#---Atmospheric model---#
# Initial parameters:
date_hour = DateTime(2023, 12, 3, 13, 30, 40) # Date and hour: (y, m, d, h, mi, s,)
initial_altitude = 700000                     # Starting altitude [m]
final_altitude = 150000                       # Final altitude [m]
#latitude= -22 * (pi / 180)                    # Latitude [rad]
#longitude = -45 * (pi / 180)                  # Longitude [rad]
F10_7a = 150                                  # F10.7: 81-day average solar flux (in units of 10^-22 W/m^2/Hz)
F10_7 = 150                                   # F10.7: Daily solar flux (in units of 10^-22 W/m^2/Hz)
Ap = 4                                        # Ap magnetic index (can be a number or a vector for daily values)


#---Latitude and longitude---#
e = 0.0001               # Eccentricity
tolerance = 10^-12       # High precision tolerance used in computing the latitude