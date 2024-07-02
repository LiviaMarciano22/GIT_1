include("functions.jl")

"""
A function computing the position vector starting from the x_mee:
    - `x::Vector`: Orbital State Vector in mee

"""

function mee2positionvector(x) # LIV
    
    s = 1 + x[4]^2 + x[5]^2                   # It is 's^2' in the pdf pg 3
    q = 1 + x[2]*cos(x[6]) + x[3]*sin(x[6])   # it is 'w' in the provided pdf, page 3 
    r = x[1]/q
    a_2 = x[4]^2 - x[5]^2                     # It is 'alpha^2' in the pdf pg 3
    pos_vec_ECI = zeros(3)                    # Position vector in ECI
    pos_vec_ECI[1] = (r/s)*(cos(x[6]) + a_2*cos(x[6]) + 2*x[4]*x[5]*sin(x[6]))  # First component of position vector
    pos_vec_ECI[2] = (r/s)*(sin(x[6]) - a_2*sin(x[6]) + 2*x[4]*x[5]*cos(x[6]))  # Second component of position vector
    pos_vec_ECI[3] = (2*r/s)*(x[4]*sin(x[6]) - x[5]+cos(x[6]))                  # Third component of position vector
    
    return pos_vec_ECI

end 

"""
A function computing Greenwich Mean Sidereal Time in degrees
"""

function calculate_gmst(date_hour) #LIV

    jd = AstroLib.juldate(date_hour)                                                                     # Calculate the Julian Day using AstroLib.jl
    t = (jd - 2451545.0) / 36525.0
    gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t^2 - t^3 / 38710000.0
    gmst = mod(gmst, 360.0)                                                                               # Normalize the result to a range from 0 to 360 degrees
    gmst_rad = gmst * (pi / 180)                                                                          # Convert to radians
   
    return gmst_rad
end 

"""
A function converting the position vector from ECI to ECEF:
      - `x::Vector`: Orbital State Vector in mee
"""

function ECI2ECEF(x, date_hour) #LIV 
    
    pos_vec_ECEF = zeros(3)
    pos_vec_ECI = mee2positionvector(x)
    gmst_rad = calculate_gmst(date_hour)

    cos_theta = cos(gmst_rad)
    sin_theta = sin(gmst_rad)

    # Rotation matrix for Z-axis rotation --> clockwise rotation
    rotation_matrix = [cos_theta sin_theta 0;
                       -sin_theta  cos_theta 0;
                       0          0         1] 

    # Perform the rotation
    pos_vec_ECEF = rotation_matrix * pos_vec_ECI  # It is doing that for each component of pos_vec_ECI
   
   return pos_vec_ECEF
   
end


"""
A function that obtains longitude and latitude from a position vector expressed in ECEF coordinates, 
as described in Algorithm 12 of 'Fundamentals of Astrodynamics and Applications' by David A. Vallado.
Thus, the method is using spherical trigonometry:
  - `x::Vector`: Orbital State Vector in mee
"""

function ECEF2LatLon(x, date_hour) #LIV
   
    pos_vec_ECI = mee2positionvector(x)
    gmst_rad = calculate_gmst(date_hour)
    pos_vec_ECEF = ECI2ECEF(pos_vec_ECI, gmst_rad)

    eq_projection_pos_sat = sqrt(pos_vec_ECEF[1]^2 + pos_vec_ECEF[2]^2)              # [km] Equatorial projection of the satellite's position vector (rdeltasat)
    longitude = asin(pos_vec_ECEF[2]/eq_projection_pos_sat)                          # [deg] Longitude --> because the vector is Earth fixed, the right ascension is also equal to the longitude 
    # To find the geodetic latitude it is needed iteration 
    distance_E_S_z = pos_vec_ECEF[3]                                                 # [km] Distance from the center of the Earth to the satellite along the z axis (rksat)
    delta = atan(distance_E_S_z/eq_projection_pos_sat)                               # [deg] First guess for the latitude 
    phi_old = delta                                                                  # [deg] Latitude (phiold)--> first guess 
    while true                                                                      
    C_m = Re/ (sqrt(1-(e^2)*(sin(phi_old))^2))                                       # [km] Radius of curvature in the meridian --> Auxiliary quantity obtained from geometrical properties of an ellipse 
    phi_new = atan((distance_E_S_z + C_m*e^2*sin(phi_old))/eq_projection_pos_sat)
    # Check if the change in latitude is within the tolerance
     if abs(phi_new - phi_old) < tolerance
        break
    end
    # Update the old latitude value
    phi_old = phi_new
    end
    latitude = phi_old                                                               # [deg] Geodetic Latitude 
    
    return latitude,longitude
end


"""
Computes the Gauss Variational Equations of motion of the Spacecraft in MEE coordinates in RTN frame
# ARGUMENTS
- `x::Vector`: Orbital State Vector in mee
- `t::Vector`: Time vector
- `u::Vector`: Control's Array
- `J2_per::Bool`: if `True` J2 perturbation active
- `Drag_per::Bool`: if `True` Drag perturbation active
"""
function dynamics(x,date_hour,t,u,J2_per,Drag_per) # LIV: aggiunto date_hour e Drag_per

   # x_mee = coe2mee(x[1:7])                    # LIV
   
    #---Auxiliar expressions---#
    A = x[4]*sin(x[6]) - x[5]*cos(x[6])
    B = x[4]*cos(x[6]) + x[5]*sin(x[6])
    s = 1 + x[4]^2 + x[5]^2
    q = 1 + x[2]*cos(x[6]) + x[3]*sin(x[6])
    r = x[1]/q
    #alpha = 1- x[4]^2 - x[5]^2
    C = 1- x[4]^2 - x[5]^2                    # it is not something defined in the pdf, it is defined to simplied the writing of the equations --> LIV
    v_r = (sqrt(μ/x[1]))*(x[2]*sin(x[6])-x[3]*cos(x[6])) #LIV --> component of velocity appearing inside the equation of the drag (page 6 file 'equazioni de implementare)
    v_t = sqrt(μ/x[1])*(1+x[2]*cos(x[6])+x[3]*sin(x[6])) #LIV --> component of velocity appearing inside the equation of the drag (page 6 file 'equazioni de implementare)
    v = sqrt(v_r^2+v_t^2)                     # Velocity magnitude --> LIV 


  # Basic gravitational perturbation forces --> LIV
   Δp_r = (u[1]*(10^-3)/x[7])*cos(u[3])*sin(u[2])
   Δp_t = (u[1]*(10^-3)/x[7])*cos(u[3])*cos(u[2])
   Δp_h = (u[1]*(10^-3)/x[7])*sin(u[3])
   
   #---Perturbations---# --> LIV
   if J2_per
       Δp_r += ((-3*μ*J2*Re^2)/(2*(r)^4))*(1-((12*(A)^2)/(s)^2))
       Δp_t += ((-12*μ*J2*Re^2)/((r)^4))*(A*B/(s)^2)
       Δp_h += ((-6*μ*J2*Re^2)/((r)^4))*(C*A/(s)^2)
   end
   
   
   if Drag_per

     x_coe = mee2coe(x[1:7])
     current_altitude = x_coe[1]-6371 
     #x_mee = coe2mee(x[1:7])                    # LIV           
     pos_vec_ECI = mee2positionvector(x)
     gmst_rad = calculate_gmst(date_hour)
     pos_vec_ECEF = ECI2ECEF(pos_vec_ECI, gmst_rad)
     lat,long = ECEF2LatLon(pos_vec_ECEF, gmst_rad)


        Atm_mod = AtmosphericModels.nrlmsise00(date_hour, current_altitude, lat, long, F10_7a, F10_7, Ap)
        
        # Update Δp_r and Δp_t with the current values
        Δp_r += ((-1/2 * Atm_mod.total_density * Areas[1] * C_d * v * v_r))
        Δp_t += ((-1/2 * Atm_mod.total_density * Areas[1] * C_d * v * v_t))
        # Δp_h remains the same for Drag
        
   end 


    #---EDO--#
    δp = sqrt(x[1]/μ)*(2*x[1]/q)*Δp_t
    δf = sqrt(x[1]/μ)*(sin(x[6])*Δp_r + (1/q)*((1+q)*cos(x[6]) + x[2])*Δp_t - (x[3]/q)*(x[4]*sin(x[6]) - x[5]*cos(x[6]))*Δp_h)
    δg = sqrt(x[1]/μ)*(-cos(x[6])*Δp_r + (1/q)*((1+q)*sin(x[6]) + x[3])*Δp_t + (x[2]/q)*(x[4]*sin(x[6]) - x[5]*cos(x[6]))*Δp_h)
    δh = sqrt(x[1]/μ)*(s*cos(x[6])/(2*q))*Δp_h
    δk = sqrt(x[1]/μ)*(s*sin(x[6])/(2*q))*Δp_h
    δL = sqrt(x[1]/μ)*(((1/q)*(x[4]*sin(x[6]) - x[5]*cos(x[6])))*Δp_h) + sqrt(x[1]*μ)*(q/x[1])^2
    δm = -u[1]/(Isp*g0)

    return [δp, δf, δg, δh, δk, δL, δm]
end
"""
Runge Kutta 4th order implementation
# ARGUMENTS
- `x::Array`: Orbital State Array
- `t::Vector`: Time Vector
- `u::Array`: Control's Array
- `f::Function`: Equations of motion
- `i::Int64`: Step index
- `J2_per::Bool`: if `True` J2 perturbation active
"""

function Rk4(x, t, u, f, i, J2_per, Drag_per) # LIV ho aggiunto drag per
    h = t[i+1] - t[i]
    k1 = f(x[i,:], t[i], u[i,:], J2_per, Drag_per)
    k2 = f(x[i,:] + k1 * h/2, t[i] + h/2, u[i,:], J2_per, Drag_per)
    k3 = f(x[i,:] + k2 * h/2, t[i] + h/2, u[i,:], J2_per, Drag_per)
    k4 = f(x[i,:] + k3 * h, t[i] + h, u[i,:], J2_per, Drag_per)
    x[i+1,:] += x[i,:] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    x_coe[i+1,:] += mee2coe(x_mee[i+1,:])
end


"""
Initialize the Arrays
# ARGUMENTS
- `x_mee::Array`: Orbital State's Array in MEE
- `x_coe::Array`: Orbital State's Array in COE
- `x0::Vector`: Initial Value's vector
- `u::Array`: Control's Array
"""
function Initialize(x_mee,x_coe, x0, u)
    x_mee[1,:] += coe2mee(x0)
    x_coe[1,:] += x0
    u[1,1] += 0
    u[1,2] = 0
    u[1,3] = 0
    Temp[1,:] .+= Temp_init
end