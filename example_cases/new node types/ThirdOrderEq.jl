"""
A node type that applies the swing equation and voltage dynmaics
to the frequency/angle dynamics.
In the following, we followed the implementation of the Extended Model/ third
order model including voltage dynamics according to K. Schmietendorf et. al.
"Self-Organized Synchronization and Voltage Stability in Networks of Synchronous
Machines", 2013.
The third order model is a special case of the fourth order model with Ed = 0
and Xd = Xd_dash.
Additionally to u, it has the internal dynamic variable ω representing the
frequencyof the rotator relative to the grid frequency Ω, i.e. the real
frequency ωr of the rotator is given as ωr=Ω+ω.
Keyword Arguments
H: inertia constant (given in [s]
P: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
D: damping coefficient, (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal Dω)
Ω: rated frequency in [1/s] of the power grid, often 2π⋅50Hz
T_d_dash: time constant of d-axis, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general explanation on time constants
X_d_dash: transient reactance of d-axis, given in [pu]
X_d: reactance of d-, given in [pu]
X_d: reactance of q-axis, given in [pu]
E_f: scaled field voltage, which, if set equal to 1.0 pu, gives 1.0 pu open-circuit terminal voltage. The physical device that provides the value of `E_f` is called the exciter (according to P. Sauer, p. 65)
"""

using PowerDynamics
using NetworkDynamics
import PowerDynamics: dimension, symbolsof, construct_vertex
import PowerDynamics: AbstractNode
import PowerDynamics: showdefinition
import Base: @__doc__

@DynamicNode ThirdOrderEq(H, P, D, Ω, E_f,T_d_dash, X_q_dash, X_d_dash, X_d) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert X_q_dash >= 0
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"
    Ω_H = Ω / (2 * H)
end [[θ,dθ],[ω, dω]] begin
    i_c = 1im * i * exp(-1im * θ)
    e_c = 1im * u * exp(-1im * θ)
    p = real(u * conj(i))
    e_d = 0                         # simpilfication from K.Schmietendorf paper
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    dθ = ω
    dω = (P - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
    de_q = (1 / T_d_dash) * (E_f - e_q + i_d * (X_d - X_d_dash))
    # -> u = e_q * exp(1im * θ)
    du = de_q * exp(1im * θ) + u * 1im * ω

    #de_c = 1im*de_q
    #du = -1im * de_c * exp(1im * θ)+ u * 1im * ω
    #dω = (P - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
end

export ThirdOrderEq
