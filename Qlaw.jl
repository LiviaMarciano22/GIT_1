
# NO CHANGES BY LIV 

"""
Computes Q-law thrust angles.
# Inputs
- `current`: osculating classical orbital elements;
- `target`: target classical orbital elements;
- `weights`: Q-law elements weights list;
- `override_ν`: true anomaly override flag;
- `ν_imposed`: imposed true anomaly (if overridden) [rad];
# Returns
- `α`: in-plane thrust angle [rad];
- `β`: out-of-plane thrust angle [rad];
- `Qdotmin`: minimum Q rate;
"""
function qlaw_angles(current, target, weights) 
    mq, nq, rq = (3, 4, 2)
    # if averaged
    #     (; a, e, i, Ω, ω, ν) = J2AveragedElements(current)
    # else

    a = current[1]
    e = current[2]
    i = current[3]
    ω = current[4]
    Ω = current[5]
    ν = current[6]

    #(a, e, i, Ω, ω, ν) .= current
    aₜ = target[1]
    eₜ = target[2]
    iₜ = target[3]
    ωₜ = target[4]
    Ωₜ = target[5]
    #(aₜ, eₜ, iₜ, Ωₜ, ωₜ) .= target

    #Wa, We, Wi, WΩ, Wω = weights

    Wa = weights[1]
    We = weights[2]
    Wi = weights[3]
    Wω = weights[4]
    WΩ = weights[5]
    km = 1000
    s = 1

    p = a * (1 - e^2)
    #println(p)
    h = √(p * μ)
    r = p / (1 + e * cos(ν))
    # Derivatives of Qa ____________________________________________________________________________________________
    if !iszero(Wa)
        dQa_da = Wa * nq * ((a - aₜ) / (aₜ * mq))^nq * (a^3 * (e + 1) / (μ * (1 - e)))^(-1.0) * (a - aₜ) *
                 (((a - aₜ) / (aₜ * mq))^nq + 1)^(1 / rq) / (4 * rq * (((a - aₜ) / (aₜ * mq))^nq + 1)) +
                 Wa * (a^3 * (e + 1) / (μ * (1 - e)))^(-1.0) * (2 * a - 2 * aₜ) * (((a - aₜ) / (aₜ * mq))^nq + 1)^(1 / rq) / 4 -
                 0.75 * Wa * (a^3 * (e + 1) / (μ * (1 - e)))^(-1.0) * (a - aₜ)^2 * (((a - aₜ) / (aₜ * mq))^nq + 1)^(1 / rq) / a
        dQa_de = Wa * μ * (a^3 * (e + 1) / (μ * (1 - e)))^(-1.0) * (1 - e) * (a - aₜ)^2 *
                 (-1.0 * a^3 / (μ * (1 - e)) - 1.0 * a^3 * (e + 1) / (μ * (1 - e)^2)) *
                 (((a - aₜ) / (aₜ * mq))^nq + 1)^(1 / rq) / (4 * a^3 * (e + 1))
        dQa_di = 0.0km^2 / s^2
        dQa_dΩ = 0.0km^2 / s^2
        dQa_dω = 0.0km^2 / s^2
    else
        dQa_da = 0.0km / s^2
        dQa_de = dQa_di = dQa_dΩ = dQa_dω = 0.0km^2 / s^2
    end
    # Derivatives of Qe ____________________________________________________________________________________________
    if !iszero(We)
        dQe_da = -0.25 * We * (a * μ * (1 - e^2)) * (e - eₜ)^2 / (a^3 * (1 - e^2)^2)
        dQe_de = 0.5 * We * e * (a * μ * (1 - e^2)) * (e - eₜ)^2 / (a^2 * (1 - e^2)^3) +
                 We * (a * μ * (1 - e^2)) * (2 * e - 2 * eₜ) / (4 * a^2 * (1 - e^2)^2)
        dQe_di = 0.0km^2 / s^2
        dQe_dΩ = 0.0km^2 / s^2
        dQe_dω = 0.0km^2 / s^2
    else
        dQe_da = 0.0km / s^2
        dQe_de = dQe_di = dQe_dΩ = dQe_dω = 0.0km^2 / s^2
    end
    # Derivatives of Qi ____________________________________________________________________________________________
    if !iszero(Wi)
        dQi_da = -Wi * (a * μ * (1 - e^2)) * (i - iₜ)^2 * (-e * abs(cos(ω)) + √(-e^2 * sin(ω)^2 + 1))^2 / (
            a^3 * (1 - e^2)^2)
        dQi_de = 2.0 * Wi * e * (a * μ * (1 - e^2)) * (i - iₜ)^2 * (-e * abs(cos(ω)) + √(-e^2 * sin(
                     ω)^2 + 1))^2 / (a^2 * (1 - e^2)^3) + Wi * (a * μ * (1 - e^2)) * (i - iₜ)^2 * (
                                                         -e * abs(cos(ω)) + √(-e^2 * sin(ω)^2 + 1)) * (-2.0 * e * (
                                                              -e^2 * sin(ω)^2 + 1)^(-0.5) * sin(ω)^2 - 2 * abs(cos(ω))) / (a^2 * (1 - e^2)^2)
        dQi_di = Wi * (a * μ * (1 - e^2)) * (2 * i - 2 * iₜ) * (-e * abs(cos(ω)) + √(-e^2 * sin(ω)^2 + 1))^2 / (
            a^2 * (1 - e^2)^2)
        dQi_dΩ = 0.0km^2 / s^2
        dQi_dω = 0.0km^2 / s^2  # TODO: complete equation
    else
        dQi_da = 0.0km / s^2
        dQi_de = dQi_di = dQi_dΩ = dQi_dω = 0.0km^2 / s^2
    end
    # Derivatives of QΩ ____________________________________________________________________________________________
    if !iszero(WΩ)
        dQΩ_da = WΩ * ((μ * acos(cos(Ω - Ωₜ))^2 * sin(i)^2 * ((1 - e^2 * cos(ω)^2)^(1 / 2) - e *
                                                                                             abs(sin(ω)))^2) / (a^2 * (e^2 - 1)))
        dQΩ_de = WΩ * ((2 * e * μ * acos(cos(Ω - Ωₜ))^2 * sin(i)^2 * ((1 - e^2 * cos(ω)^2)^(1 / 2) - e
                                                                                                     *
                                                                                                     abs(sin(ω)))^2) / (a * (e^2 - 1)^2) + (2 * μ * acos(cos(Ω - Ωₜ))^2 * sin(i)^2 * ((1
                                                                                                                                                                                       -
                                                                                                                                                                                       e^2 * cos(ω)^2)^(1 / 2) - e * abs(sin(ω))) * (abs(sin(ω)) + (e * cos(ω)^2) / (1 - e^2 * cos(ω)^2)^(
                                                                                                                                                1 / 2))) / (a * (e^2 - 1)))
        dQΩ_di = -WΩ * (2 * μ * acos(cos(Ω - Ωₜ))^2 * cos(i) * sin(i) * ((1 - e^2 * cos(ω)^2)^(1 / 2) - e *
                                                                                                        abs(sin(ω)))^2) / (a * (e^2 - 1))
        dQΩ_dΩ = -WΩ * (2 * μ * sin(Ω - Ωₜ) * acos(cos(Ω - Ωₜ)) * sin(i)^2 * ((1 - e^2 * cos(ω)^2)^(1 / 2) -
                                                                              e * abs(sin(ω)))^2) / (a * (e^2 - 1) * (1 - cos(Ω - Ωₜ)^2)^(1 / 2))
        dQΩ_dω = 0.0  # TODO: complete equation
    else
        dQΩ_da = 0.0km / s^2
        dQΩ_de = dQΩ_di = dQΩ_dΩ = dQΩ_dω = 0.0km^2 / s^2
    end
    if !iszero(Wω)
        dQω_da = 0.0km / s^2
        dQω_de = dQω_di = dQω_dΩ = dQω_dω = 0.0km^2 / s^2 # TODO: complete equations
    else
        dQω_da = 0.0km / s^2
        dQω_de = dQω_di = dQω_dΩ = dQω_dω = 0.0km^2 / s^2
    end
    # Full derivatives per element _________________________________________________________________________________
    dQ_da = dQa_da
    dQ_da += dQe_da
    dQ_da += dQi_da
    dQ_da += dQΩ_da
    dQ_da += dQω_da
    dQ_de = dQa_de
    dQ_de += dQe_de
    dQ_de += dQi_de + dQΩ_de + dQω_de
    dQ_di = dQa_di + dQe_di + dQi_di + dQΩ_di + dQω_di
    dQ_dΩ = dQa_dΩ + dQe_dΩ + dQi_dΩ + dQΩ_dΩ + dQω_dΩ
    dQ_dω = dQa_dω + dQe_dω + dQi_dω + dQΩ_dω + dQω_dω
    # Variation with respect to thrust _____________________________________________________________________________
    d_da_dFt = (2 * (a^2) / h) * (1 + e * cos(ν))
    d_da_dFr = (2 * (a^2) / h) * e * sin(ν)
    d_da_dFn = 0.0s
    d_de_dFt = 1 / h * ((p + r) * cos(ν) + r * e)
    d_de_dFr = 1 / h * p * sin(ν)
    d_de_dFn = 0.0s / km
    d_di_dFt = 0.0s / km
    d_di_dFr = 0.0s / km
    d_di_dFn = r * cos(ν + ω) / h
    d_dΩ_dFt = 0.0s / km
    d_dΩ_dFr = 0.0s / km
    d_dΩ_dFn = r * sin(ν + ω) / (h * sin(i))
    d_dω_dFt = 0.0s / km # TODO: complete equation?
    d_dω_dFr = 0.0s / km  # TODO: complete equation?
    d_dω_dFn = 0.0s / km # TODO: complete equation?
    # D1, D2 and D3 ________________________________________________________________________________________________
    D1 = dQ_da * d_da_dFt + dQ_de * d_de_dFt + dQ_di * d_di_dFt + dQ_dΩ * d_dΩ_dFt + dQ_dω * d_dω_dFt
    D2 = dQ_da * d_da_dFr + dQ_de * d_de_dFr + dQ_di * d_di_dFr + dQ_dΩ * d_dΩ_dFr + dQ_dω * d_dω_dFr
    D3 = dQ_da * d_da_dFn + dQ_de * d_de_dFn + dQ_di * d_di_dFn + dQ_dΩ * d_dΩ_dFn + dQ_dω * d_dω_dFn
    # Adding Qdotmin here, means that is not necessary to recompute it in another function ________________________
    Qdotmin = -√(D1^2 + D2^2 + D3^2)
    α = atan(-D2, -D1) # -α_regularization, α_regularization)  # Computation and storage of α thrust angle
    β = atan((-D3) / (√(D1^2 + D2^2)))  # Computation and storage of β angle
    # @debug "$D1 $D2 $D3 $α $β $Qdotmin"
    return (α, β)
end