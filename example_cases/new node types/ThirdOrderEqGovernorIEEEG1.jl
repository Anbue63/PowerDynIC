"""
A node type that applies the swing equation and voltage dynmaics
to the frequency/angle dynamics.
In the following, we followed the implementation of the Ex11tended Model/ third
order model including voltage dynamics according to K. Schmietendorf et. al.
"Self-Organized Synchronization and Voltage Stability in Networks of Synchronous
Machines", 2013.
The third order model is a special case of the fourth order model with Ed = 0
and x11d = x11d_dash.

In addition to the third order model the  Type IEEEG1 governor was implemented.
It is used for tandem compound, double reheat steam turbine systems.
The type was implemented according to "Dynamic Models for Turbine-Governors in
Power System Studies", IEEE PES, 2013, Page 2 - 2

Additionally to u, it has the internal dynamic variable ω representing the
frequencyof the rotator relative to the grid frequency Ω, i.e. the real
frequency ωr of the rotator is given as ωr=Ω+ω.

Keyword Arguments

Third Order Model
- H: inertia constant (given in [s]
- P: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- D: damping coefficient, (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal Dω)
- Ω: rated frequency in [1/s] of the power grid, often 2π⋅50Hz
- T_d_dash: time constant of d-ax11is, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general ex11planation on time constants
- X_d_dash: transient reactance of d-axis, given in [pu]
- X_q_dash: transient reactance of q-axis, given in [pu]
- X_d: reactance of d-, given in [pu]
- X_d: reactance of q-ax11is, given in [pu]
- E_f: scaled field voltage, which, if set equal to 1.0 pu, gives 1.0 pu open-circuit terminal voltage. The physical device that provides the value of `E_f` is called the ex11citer (according to P. Sauer, p. 65)

Governer Type
- P0:
- Pmax: Max power limit imposed by valve or gate control [pu]
- Pmin: Min power limit imposed by valve or gate control [pu]
- Pup: Max main control valve rate of change [pu/s]
- Pdown: Min main control valve rate of change [pu/s]
- T1: Controller lag compensation [s]
- T2: Controller lead compensation [s]
- T3: Governor lag [s]
- K: Total effective speed-governing system gain [pu]
"""

using PowerDynamics
using NetworkDynamics
import PowerDynamics: dimension, symbolsof, construct_vertex
import PowerDynamics: AbstractNode
import Base: @__doc__

@DynamicNode ThirdOrderEqGovernorIEEEGone(H, D, Ω, E_f,T_d_dash, X_d_dash, X_q_dash, X_d, P0, Pmax, Pmin, Pup, Pdown, T1, T2, T3, K) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-ax11is (x11_d_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-ax11is (x11_d_dash) should be >=0"
    Ω_H = Ω / (2 * H)
end [[θ,dθ],[ω, dω], [Pm, dPm], [x1, dx1], [z, dz], [P, dP]] begin
    i_c = 1im * i * exp(-1im * θ)
    e_c = 1im * u * exp(-1im * θ)
    p = real(u * conj(i))
    e_d = 0                         # simpilfication from K.Schmietendorf paper
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    # Governor Type IEEEG1
    dx1 = K * (-1 / T1 * x1 + (1 - T2 /  T1) * ω) # Block Input

    dP = (1 / T1) * x1 + (T2 / T1) * ω

    y = (1 / T3) * (P0 - P - Pm)                  # Block Output
    y_temp = y                                    # temorary variable

    # Limiting the valve rate of change
    if y > Pup
        y_temp = Pup
    end

    if y < Pdown
        y_temp = Pdown
    end

    dz = y_temp
    dPm = y_temp

    # Limiting the power imposed by the valve or gate control
    if z > Pmax
        dPm = (1 - Pmax) * dPm
    end

    if z < Pmin
        dPm = (1 - Pmin) * dPm
    end

    # Third Order Model
    dθ = ω
    dω = (Pm - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
    de_q = (1 / T_d_dash) * (E_f - e_q + i_d * (X_d - X_d_dash))
    # -> u = e_q * exp(1im * θ)
    du = de_q * exp(1im * θ) + u * 1im * ω
end

export ThirdOrderEqGovernorIEEEGone
