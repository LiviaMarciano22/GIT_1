# NO CHANGES BY LIV 

#----Functions-----#

#NODE DEFINITION 
#   1)  [1  0 0] radio vector (r)
#   2)  [-1 0 0]
#   3)  [0  1 0] Tangential direction ((r x v) x r) (Thrust direction)
#   4)  [0 -1 0] 
#   5)  [0  0 1] Normal to the orbital plane (r x v)
#   6)  [0 0 -1]

"""
Calculates the angles between Earth, sun and surface's normals
# ARGUMENTS
- `Earth_2_normal::vector`: Contains the angles between the Spacecraft's position vector and the normals to the surfaces [radians].
- `Sun_2_normal::vector`: Contains the angles between the Sun rays and the normals to the surfaces [radians].
- `Sun_2_Earth::vector`: Contains the angles between the Sun rays and the Spacecraft's position vector [radians].
- `r_sun::vector`: Position of the sun in ECI frame [km].
- `x_sc::vector`: Spacecraft's State vector [km,km/s].
- `node::int`: Node number index
- `α::float64`: Thrust yaw angle [radians].
- `β::float64`: Thrust pitch angle [radians].
"""
function Incidence_Angles(Earth_2_normal, Sun_2_normal, Sun_2_Earth, r_sun, x_sc, node, α, β)
    #Local axis definition
    xBody = [1;0;0]
    yBody = [0;1;0] 
    zBody = [0;0;1]
    n_vectors= [xBody , -xBody , yBody , -yBody , zBody , -zBody, -yBody]

    n_RTN = (Rz(α)*Rx(-β))*n_vectors[node] #normal Vector in RTN frame
    n_ECI = RTNtoECI(x_sc)*n_RTN

    r_sc=x_sc[1:3] #Spacecraft radio vector
    S = r_sun - r_sc
    Earth_2_normal[node] = acos(round(dot(-r_sc,n_ECI)/(norm(r_sc)*norm(n_ECI)),digits=5))
    Sun_2_normal[node] = acos(round(dot(n_ECI,S)/(norm(n_ECI)*norm(S)),digits=5))
    Sun_2_Earth[node] = acos(round(dot(r_sc,r_sun)/(norm(r_sc)*norm(r_sun)),digits=5))
    return nothing
end
"""
Calculates the view factor Fe between the Earth and the Spacecraft's surfaces.
# ARGUMENTS
- `Earth_2_normal`::vector: Contains the angles between the Spacecraft's position vector and the normals to the surfaces [radians].
- `r_sc::vector`: Position of the Spacecraft in ECI frame [km].
- `node::int`: Node number index.
"""
function Earth_ViewFactor(Earth_2_normal, r_sc, node)
    ϕ = Earth_2_normal[node]
    R_sc = norm(r_sc)
    Re = 6378
    H = R_sc/Re
    Φm = asin(1/H)
    fact = sqrt(H^2 -1)

    if ϕ ≤ pi/2 - Φm
        return Fe = cos(ϕ)/(H^2)

    elseif (ϕ ≤ pi/2 + Φm)   && (ϕ > pi/2 - Φm)
        t1 = 0.5 - (1/pi)*asin(fact/(H*sin(ϕ)))
        t2 = (1/(pi*H^2))*(cos(ϕ)*acos(-fact*cot(ϕ)))
        t3 = fact*sqrt(1-(H^2)*(cos(ϕ))^2 )*(1/(pi*H^2))
        return Fe = t1 + t2 -t3
    else
        return Fe = 0
    end
end
"""
Calculates the Solar heatflux and updates the heat Array.
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `Sun_2_normal::vector`: Contains the angles between the Sun rays and the Spacecraft's position vector [radians].
- `r_sc::vector`: Position of the Spacecraft in ECI frame [km].
- `r_sun::vector`: Position of the Sun in ECI frame [km].
- `Areas::vector`: Contains the area of each surface/node [m^2].
- `Q_sun::float64`: The solar constant [W/m^2].
- `alpha::float64`: The solar absorptivity factor.
- `node::int`: Node number index.
- `j::int`: Time index.
"""
function Solar_heatflux(heat, Sun_2_normal, r_sc, r_sun, Areas, Q_sun, alpha, node, j)
    s = Step_eclipse(r_sc, r_sun)
    if Sun_2_normal[node]>= pi/2
        heat[j, node]+=0
    else
        heat[j,node]+= Q_sun*alpha[node]*Areas[node]*cos(Sun_2_normal[node])*s
    end
    return nothing
end
"""
Calculates the heatflux due to the Earth and updates heat Array.
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `Areas::vector`: Contains the area of each surface/node [m^2].
- `Fe::float64`: View factor between Earth and Spacecraft's surfaces.
- `Q_ir::float64`: The Earth IR radiation heat constant [W/m^2].
- `epsilon::float64`: Surface's IR emissivity.
- `node::int`: Node number index.
- `j::int`: Time index.
"""
function Earth_heatflux(heat, Areas, Fe,  Q_ir, epsilon, node, j)
    heat[j,node]+= Q_ir*Areas[node]*epsilon[node]*Fe
    return nothing
end
"""
Calculates the heatflux due to Albedo effect and updates the heat Array.
# ARGUMENTOS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `Sun_2_Earth::vector`: Contains the angles between the Sun rays and the Spacecraft's position vector [radians].
- `Q_sun::Float`: The solar constant [W/m^2]
- `r_sc::vector`: Position of the Spacecraft in ECI frame [km].
- `r_sun::vector`: Position of the Sun in ECI frame [km].
- `alpha::float64`: The solar absorptivity factor.
- `af::float64`: Albedo factor.
- `Fe::float64`: View factor between Earth and Spacecraft's surfaces.
- `node::int`: Node number index.
- `j::int`: Time index.
"""
function Albedo_heatflux(heat, Sun_2_Earth, Q_sun, r_sc, r_sun, alpha, af, Areas, Fe, node,j)
    θ = Sun_2_Earth[node]
    s = Step_eclipse(r_sc, r_sun)
    if θ ≥ pi/2
        heat[j,node]+= 0
    else
        heat[j,node]+= Q_sun*Areas[node]*af*alpha[node]*Fe*cos(θ)*s
    end
    return nothing
end
"""
Calculates the conductive heatflux for a given node at a given time
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `k::int`: Contact conductance coefficient [W/k].
- `Cond_area::Array`: Contains every conductive interface area for each node
- `Temp::Array`: Contains every temperature per node at every Time
- `node::int`:  Node number index
- `j::int`: Time index
"""
function Conductive_heatflux(heat, k, D, L, w, Temp, node, j)
    if node ≤ 6
        index=[1,2,3,4,5,6]
        if node % 2 == 0
            splice!(index,node-1)
        else
            splice!(index,node+1)
        end
        if node==3 || node==4
           Area_cond = [D*w, D*w, D*w, D*w, D*w, D*w]
        else
           Area_cond = [L*w, L*w, D*w, D*w, L*w, L*w]
        end
        for i in index
            heat[j,node]+=k*Area_cond[i]*(Temp[j,i]-Temp[j,node])
        end
        if node == 4
            heat[j,node]+=(Temp[j,7]-Temp[j,node])
        end
    else
        if node == 7
            heat[j,node]+=(Temp[j,4]-Temp[j,node])
            heat[j,node]+=kₜ*Aₜ*(Temp[j,8]-Temp[j,node])

        elseif node == 8
            heat[j,node]+=kₜ*Aₜ*(Temp[j,7]-Temp[j,node])
        end
    end
    return nothing
end
"""
Calculates the heat radiated by the spacecraft to the Space
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `epsilon:float64`: Surface IR emissivity
- `sigma::float64`: Surface absorptivity
- `Areas::Array`: Array containing each node superficial Area
- `Temp::Array`: Contains every temperatura per node at every Time
- `node::int`: Node number index
- `j::int`: Time index
"""
function Radiative_heatflux(heat, epsilon, sigma, Areas, Temp, node, j)
    heat[j,node]-=sigma*epsilon[node]*Areas[node]*(Temp[j,node])^4
    return nothing
end
"""
    Computes the temperature in the next time step
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `Temp::Array`: Contains every temperature per node at every Time [k].
- `C::Float`: Thermal capacity [J/K].
- `dt::float64`: Time step
- `node::int`: Node number index
- `j::int`: Time index
"""
function Temperature_update(heat, Temp,  C, dt, node, j)
    if j<length(t)
         Temp[j+1,node]+= Temp[j,node] + (dt/C[node])*heat[j,node]
    end
    return nothing
end

"""
Plots the thermal profile of every node.
"""
function Thermal_Plot(Time, Temp)

    Plot1 =Plots.plot(Time./3600, Temp[:,1], label="node 1", ylabel="Temp [°C]", xformatter=(_...) -> "")
    Plot2 =Plots.plot(Time./3600, Temp[:,2], label="node 2",  xformatter=(_...) -> "")
    Plot3 =Plots.plot(Time./3600, Temp[:,3], label="node 3", ylabel="Temp [°C]", xformatter=(_...) -> "")
    Plot4 =Plots.plot(Time./3600, Temp[:,4], label="node 4",  xformatter=(_...) -> "")
    Plot5 =Plots.plot(Time./3600, Temp[:,5], label="node 5", ylabel="Temp [°C]",  xformatter=(_...) -> "")
    Plot6 =Plots.plot(Time./3600, Temp[:,6], label="node 6",   xformatter=(_...) -> "")
    Plot7 =Plots.plot(Time./3600, Temp[:,7], label="thruster", xlabel="Time [h]",ylabel="Temp [°C]",)
    Plot8 =Plots.plot(Time./3600, Temp[:,8], label="Res+PPU", xlabel="Time [h]")
    #Plots.link_xaxes!(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6)
    temperatures = Plots.plot(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6,Plot7,Plot8, layout=Plots.grid(4,2, heights=[0.25,0.25,0.25,0.25], widths=[0.5,0.5])
    ,xrotation=0, plot_title="Temperature profile",link=:x)

    Plots.plot(t_updt,Temp[:,8],xlabel="Time [s]", ylabel="Temp [°C]", title="Temperature at RES+PPU")
    #Plots.hline!([Temp_min-273.15], color=:red, linestyle=:dash, label="Lower Limit")
    PPU_temp=Plots.hline!([Temp_max-273.15], color=:red, linestyle=:dash, label="Upper Limit")

    Plots.plot(t_updt,Temp[:,7],xlabel="Time [s]", ylabel="Temp [°C]", title="Temperature at Thruster")
    #Plots.hline!([Temp_inf-273.15], color=:red, linestyle=:dash, label="Lower Limit")
    thruster_temp=Plots.hline!([Temp_sup-273.15], color=:red, linestyle=:dash, label="Upper Limit")

    Plots.savefig(thruster_temp,"TFG_final_code/thermal_plots/thruster_temp.png")
    Plots.savefig(temperatures,"TFG_final_code/thermal_plots/temperature_profiles.png")
    Plots.savefig(PPU_temp,"TFG_final_code/thermal_plots/PPU_temp.png")

end
"""
Computes the heat in the actual step time and calculates the temperature for the next step
# ARGUMENTS
- `heat::Array`: Array containing the total heatflux for every discrete time in every node [W].
- `Temp::Array`: Contains every temperature per node at every Time [k].
- `x_coe::Vector`: Orbital State Vector in COE 
- `angles::Tuple`: Control angles in rad
- `r_sun::Vector`: Sun position vector [Km]
- `j::Int64`: Step time index
"""
function thermal_model(heat, Temp, x_coe, angles, r_sun, controls, j)
    x_sc = Perifocal2ECI(x_coe)
    r_sc = x_sc[1:3]
    α = angles[1]
    β = angles[2]
    if i==1
        Thrust = controls[i,1]
    else
        Thrust = controls[i-1,1]
    end
    for node=1:8
        if node ≤ 7
           Incidence_Angles(Earth_2_normal, Sun_2_normal, Sun_2_Earth, r_sun, x_sc, node, α, β)
           Fe = Earth_ViewFactor(Earth_2_normal, r_sc, node)
           Solar_heatflux(heat, Sun_2_normal, r_sc, r_sun, Areas, Q_sun, alpha, node, j)
           Earth_heatflux(heat, Areas, Fe,  Q_ir, epsilon, node, j)
           Albedo_heatflux(heat, Sun_2_Earth, Q_sun, r_sc, r_sun, alpha, af, Areas, Fe, node,j)
           Conductive_heatflux(heat, k, D, L, W, Temp, node, j)
           Radiative_heatflux(heat, epsilon, sigma, Areas, Temp, node, j)
           if node==7 && Thrust!=0
                heat[j,node] += (P_out_ppu - Pd)
           end  
        else
            Conductive_heatflux(heat, k, D, L, W, Temp, node, j)
            if Thrust == 0.0
                heat[j,node]-=(P_in - P_out_ppu)*0.3
            elseif  Thrust != 0
                heat[j,node]+= (P_in - P_out_ppu)
            end
        end
        Temperature_update(heat, Temp, C, Δt, node, j)
    end
end

# """
# Calculates the heat generated by the thruster_temp
# # OUTPUT
# - `Heat_sc::Float64`: Thruster's heat in W
# """
# function ThrusterHeat()
#     Ib = Id*0.7 #[A] Discharg intensity
#     Iw = 0.27*Ib #[A] Wall current
#     V_ex = Isp/g0 #[m/s] Exhaust Velocity
#     P_ib = (0.5*m_i*V_ex^2 + (5/2)*K_boltz*T_ip + q_e*V_ion)*(Ib/q_e) # Ions in the plume
#     P_eb = ((5/2)*K_boltz*T_ep + q_e*Φ)*(Ib/q_e) # Electrons in the plume
#     P_beam = P_ib + P_eb #[W] Total Beam power
#     P_rad = 2.2*q_e*V_ion*(Ib+Iw)/q_e #[W] Radiation power loss
#     P_ex = P_beam + P_rad #[W] Total exhaust power
#     Heat_sc = Pd - P_ex
#     return Heat_sc
# end