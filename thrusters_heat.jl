# NO CHANGES BY LIV 


T_ip = 2.895e4 #k Temperatura del Ion
V_ex = 600/9.81 #m/s Velocidad de salida
K_boltz = 1.380649e-23 #J/K constante de Boltzman
m_i = 2.1771e-25 #Kg masa del ion
q_e = 1.602176634e-19 #C carga del electron
V_ion = 12.1 #eV potencial de ionizacion

Ib = 0.15*0.7 #Intensidad de descarga

T_ep = 10386.94 #k Temperatura del electron
Φ = 18 #V diferencia de potencial en el catodo
Iw = 0.27*Ib
P_ib = (0.5*m_i*V_ex^2 + (5/2)*K_boltz*T_ip + q_e*V_ion)*(Ib/q_e)
P_eb = ((5/2)*K_boltz*T_ep + q_e*Φ)*(Ib/q_e)
P_beam = P_ib + P_eb
P_rad = 2.2*q_e*V_ion*(Ib+Iw)/q_e
P_ex = P_beam + P_rad
