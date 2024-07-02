include("thermal_functions.jl")

function plotting360(t, x, u)
    #example = CSV.read("TFG_final_code/res_eclipse.txt", DataFrame)
    example = readdlm("TFG_final_code/res_eclipse.txt",',',Float64)
    SMA =  Plots.plot([t ./86400, example[:,1]./86400], [x[:,1],example[:,2]], legend = false)
    ECC =  Plots.plot([t ./86400, example[:,1]./86400], [x[:,2],example[:,3]], legend = false)
    INC =  Plots.plot([t ./86400, example[:,1]./86400], [rad2deg.(mod2pi.(x[:,3])),rad2deg.(mod2pi.(example[:,4]))], legend = false) 
    AOP =  Plots.plot([t ./86400, example[:,1]./86400], [rad2deg.(mod2pi.(x[:,4])),rad2deg.(mod2pi.(example[:,6]))], legend = false)
    RAAN = Plots.plot([t ./86400, example[:,1]./86400], [rad2deg.(mod2pi.(x[:,5])),rad2deg.(mod2pi.(example[:,5]))], legend = false)
    TA   = Plots.plot([t ./86400, example[:,1]./86400], [rad2deg.(mod2pi.(x[:,6])),rad2deg.(mod2pi.(example[:,7]))], legend = false)
    mass = Plots.plot([t ./86400, example[:,1]./86400], [x[:,7], example[:,9]], legend = false)
    Thrust = Plots.plot([t./3600, example[:,1]./3600], [u[:,1], example[:,8]], legend = false, ylim=[0.001,0.004])
    alpha =  Plots.plot([t./3600, example[:,1]./3600], [rad2deg.((u[:,2])),rad2deg.((example[:,10]))], legend = false)
    beta =   Plots.plot([t./3600, example[:,1]./3600], [rad2deg.((u[:,3])),rad2deg.((example[:,11]))], legend = false)
    COE = Plots.plot(SMA, ECC, INC, AOP, RAAN, TA, layout = Plots.grid(3, 2, heights=[0.33 ,0.33, 0.33], widths=[0.5,0.5]), ylabel=["SMA (Km)" "ECC (-)" "INC (deg)" "AOP (deg)" "RAAN (deg)" "TA (deg)"],xlabel="Time (Days)",title="Classical Orbital Elements")
    cont = Plots.plot(Thrust, mass, alpha, beta, layout = Plots.grid(2,2, heights=[0.5,0.5], widths=[0.5,0.5]), ylabel=["Thrust (N)" "Mass (Kg)" "In-plane (deg)" "Out-of-plane (deg)"],xlabel="Time (hours)", title="Controls")
    Plots.savefig(COE,"TFG_final_code/transfer_plots/COE_eclipse.png")
    Plots.savefig(cont,"TFG_final_code/transfer_plots/controls_eclipse.png")
end

function plotting(t, x, u)

    SMA =  Plots.plot([t ./86400], [x[:,1]], legend = false, xformatter=(_...) -> "")
    ECC =  Plots.plot([t ./86400], [x[:,2]], legend = false, xformatter=(_...) -> "")
    INC =  Plots.plot([t ./86400], [rad2deg.(mod2pi.(x[:,3]))], legend = false, xformatter=(_...) -> "") 
    AOP =  Plots.plot([t ./86400], [rad2deg.(mod2pi.(x[:,4]))], legend = false, xformatter=(_...) -> "")
    RAAN = Plots.plot([t ./86400], [rad2deg.(mod2pi.(x[:,5]))], legend = false, xlabel="Time (Days)")
    TA   = Plots.plot([t ./86400], [rad2deg.(mod2pi.(x[:,6]))], legend = false, xlabel="Time (Days)")
    mass = Plots.plot([t ./3600], [x[:,7]], legend = false, xformatter=(_...) -> "")
    Thrust = Plots.plot([t./3600], [u[:,1]], legend = false, ylim=[0.0,0.007], xformatter=(_...) -> "")
    alpha =  Plots.plot([t./3600], [rad2deg.((u[:,2]))], legend = false,xlabel="Time (hours)")
    beta =   Plots.plot([t./3600], [rad2deg.((u[:,3]))], legend = false,xlabel="Time (hours)")
    COE = Plots.plot(SMA, ECC, INC, AOP, RAAN, TA, layout = Plots.grid(3, 2, heights=[0.33 ,0.33, 0.33], widths=[0.5,0.5]), ylabel=["SMA (Km)" "ECC (-)" "INC (deg)" "AOP (deg)" "RAAN (deg)" "TA (deg)"],plot_title="Classical Orbital Elements",link=:x)
    cont = Plots.plot(Thrust, mass, alpha, beta, layout = Plots.grid(2,2, heights=[0.5,0.5], widths=[0.5,0.5]), ylabel=["Thrust (N)" "Mass (Kg)" "In-plane (deg)" "Out-of-plane (deg)"], plot_title="Controls",link=:x)

    Plots.savefig(COE,"TFG_final_code/thermal_plots/COE_thermal.png")
    Plots.savefig(cont,"TFG_final_code/thermal_plots/controls_thermal.png")
end
"""
Transform Modified Equinoctial Elements to Clasical Orbital Elements
# ARGUMENTS
- `x_mee::Vector`: Orbital state vector in MEE
# RETURN
- `x_coe::Vector`: Orbital state vector in COE
"""
function mee2coe(x_mee)
    x_coe = zeros(length(x_mee))
    x_coe[1] = (x_mee[1])/(1 - (x_mee[2])^2 - (x_mee[3])^2) #SMA
    x_coe[2] = sqrt(x_mee[2]^2 + x_mee[3]^2) #ECC
    x_coe[3] = 2*atan(sqrt(x_mee[4]^2 + x_mee[5]^2)) #INC
    x_coe[4] = atan(x_mee[3],x_mee[2]) - atan(x_mee[5],x_mee[4]) #Argument periapsis
    x_coe[5] = atan(x_mee[5], x_mee[4]) #RAAN
    x_coe[6] = x_mee[6]-atan(x_mee[3],x_mee[2]) #True Anomaly
    x_coe[7] = x_mee[7]
    return x_coe
end
"""
Transforms Clasical Orbital Elements to Modified Equinoctial Elements
# ARGUMENTS
- `x_coe::Vector`: Orbital state vector in COE
# RETURN
- `x_mee::Vector`: Orbital state vector in MEE
"""
function coe2mee(x_coe)
    x_mee = zeros(length(x_coe))
    x_mee[1] = x_coe[1]*(1-x_coe[2]^2)
    x_mee[2] = x_coe[2]*cos(x_coe[4] + x_coe[5])
    x_mee[3] = x_coe[2]*sin(x_coe[4] + x_coe[5])
    x_mee[4] = tan(x_coe[3]/2)*cos(x_coe[5])
    x_mee[5] = tan(x_coe[3]/2)*sin(x_coe[5])
    x_mee[6] = x_coe[5] + x_coe[4] + x_coe[6]
    x_mee[7] = x_coe[7]
    return x_mee
end

"""
Updates the control's Array in each step of time
# ARGUMENTS
- `x::Vector`: Orbital state vector in coe
- `controls::Array`: Control's Array
- `angles::Tuple`: Control angles in rad
- `i::Int64`: Time index
"""
function update_control(x, controls,angles,Temp, r_sun, EclipseC, i)
    x_sc = Perifocal2ECI(x)
    if EclipseC
        if Step_eclipse(x_sc[1:3],r_sun) == 1 
            if Temp[i,4] > Temp_max
               controls[i,1] += 0
            else
               controls[i,1] += Tmax
            end
            if (Temp[i,4] < Temp_min) && (controls[i,1]==0.0)
                controls[i,1] += Tmax
            end
        else
            if Temp[i,7] < Temp_min
                controls[i,1] += Tmax
            end
        end
    else
        if Temp_min ≤ Temp[i,8] ≤ Temp_max

            controls[i,1] += Tmax

        elseif Temp[i,8] ≤ Temp_min

            controls[i,1] += Tmax

        end
    end
    controls[i,2] = angles[1]
    controls[i,3] = angles[2]
end
"""
Adjust the Array dimension to the final maneuver time
# ARGUMENTS
- `t::Vector`: Time vector
- `x::Array`: Orbital State Array
- `u::Array`: Control's Array
- `t::Vector`: End time index
# RETURNS
- `x_updt::Array`: Adjusted Orbital State Array
- `u_updt::Array`: Adjusted Control Array
- `t_updt::Vector`: Adjusted Time Vector
"""
function updt_size(t,x,u,i)
    x_updt = x[1:i,:] 
    t_updt = t[1:i] 
    u_updt = u[1:i,:] 
    return x_updt, t_updt, u_updt
end

"""
Calculates the spacecraft's position and velocity in Perifocal frame and returns them in ECI frame
# ARGUMENTS
- `x::Array`: Orbital state vector 
# RETURNS
- `x_sc::Vector`: Spacecraft's state vector [Km , Km/s]
"""
function Perifocal2ECI(x)
    p = x[1]*(1-x[2]^2)
    r_perifocal = (p/(1+x[2]*cos(x[6]))*[cos(x[6]); sin(x[6]); 0])
    v_perifocal = (μ/sqrt(p*μ))*[-sin(x[6]); x[2]+cos(x[6]); 0]
    Rot_mat = Transpose(Rz(x[4])*Rx(x[3])*Rz(x[5]))
    r_sc = Rot_mat*r_perifocal
    v_sc = Rot_mat*v_perifocal
    x_sc = vcat(r_sc,v_sc)
    return x_sc
end

"""
Calculates whether the spacecraft is in eclipse, returns 1 if sunlight and 0 if eclipse.
# ARGUMENTS
- `r_sc::vector`: Position of the Spacecraft in ECI frame [Km].
- `r_sun::vector`: Position of the Sun in ECI frame [Unitary].
# RETURNS
- `0`: If in eclipse
- `1`: If in sunlight
"""
function Step_eclipse(r_sc, r_sun)
    r_par = dot(r_sun,r_sc)/norm(r_sun)
    r= norm(r_sc)
    if r_par + sqrt(r^2 - Re^2) ≤ 0
        return 0
    else
        return 1

    end
end
"""
Rotation matrix in X axis
"""
function Rx(θ)
    Rx = [1    0     0;
          0 cos(θ)  sin(θ);
          0 -sin(θ) cos(θ)]
    return Rx
end
"""
Rotation matrix in Y axis
"""
function Ry(θ)
    Ry = [cos(θ)    0    -sin(θ);
          0         1      0;
          sin(θ)    0    cos(θ)]
    return Ry
end
"""
Rotation matrix in Z axis
"""
function Rz(θ)
    Rz = [cos(θ)  sin(θ)  0;
          -sin(θ) cos(θ)  0;
          0         0     1]
    return Rz
end
"""
    Computes the rotation matrix from RTN frame to ECI
# ARGUMENTS
- `x_sc::vector`: Spacecraft's State Vector [Km, Km/s]
# RETURNS
- `Q::Matrix`: Rotation matrix Q RTN -> ECI
"""
function RTNtoECI(x_sc)
    r_sc = x_sc[1:3]
    v_sc = x_sc[4:6]

    i_r = r_sc./norm(r_sc)
    i_h = (cross(r_sc,v_sc))/(norm(r_sc)*norm(v_sc))
    i_t = (cross(i_h,r_sc))/(norm(i_h)*norm(r_sc))
    Q = hcat(i_r, i_t, i_h)
    return Q
end 
"""
Plots the orbit in 3D
"""
function Plot3d_Orbit(x_updt)
    position = zeros((i,3))
    for j=1:i
       x_sc=Perifocal2ECI(x_updt[j,:])
       position[j,1] = x_sc[1]
       position[j,2] = x_sc[2]
       position[j,3] = x_sc[3]
    end
Plots.plot3d(position[:,1],position[:,2],position[:,3], xlabel="X (Km)", ylabel="Y (Km)", zlabel="Z (Km)",title="Final Orbit", legend=false)
end
"""
Computes the actual time in the maneuver and updates the sun position vector
# ARGUMENTS
- `real_time::Float64`: Time counter in seconds
# RETURN
- `r_sun::Vector`: Sun position vector at actual Time
- `real_time::Float64`: Actual time
"""
function SunPositionUpdate(real_time)
    real_time += Δt
    Day = Initial_Day
    if real_time ≥ 86400
        real_time = 0
        Day += 1
    end
    #---Sun position--#
    jd= juldate(DateTime(Initial_Year, Initial_Month, Day)) #year/month/day
    vector = xyz(jd,2000)
    r_sun=[vector[1]; vector[2]; vector[3]]
    return r_sun, real_time
end