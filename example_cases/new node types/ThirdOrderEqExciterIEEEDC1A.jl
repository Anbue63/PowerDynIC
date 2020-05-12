"""
A node type that applies the swing equation and voltage dynmaics
to the frequency/angle dynamics.
In the following, we followed the implementation of the Extended Model/ third
order model including voltage dynamics according to K. Schmietendorf et. al.
"Self-Organized Synchronization and Voltage Stability in Networks of Synchronous
Machines", 2013.
The third order model is a special case of the fourth order model with e_d = 0
and X_d = X_d_dash.

In addition to the third order model the Type IEEEDC1A exciter was implemented.

Additionally to u, it has the internal dynamic variable ω representing the
frequencyof the rotator relative to the grid frequency Ω, i.e. the real
frequency ωr of the rotator is given as ωr=Ω+ω.

Keyword Arguments

Third Order Model
- H: inertia constant (given in [s]
- P: active (real) power output, also called the mechanical torque applied to the shaft, given in [pu]
- D: damping coefficient, (given in [s], see P. Sauer, eq. (5.156) where the damping torque is equal Dω)
- Ω: rated frequency in [1/s] of the power grid, often 2π⋅50Hz
- T_d_dash: time constant of d-axis, given in [s], see P. Sauer, chapter 3.7, p. 54 for a general explanation on time constants
- X_d_dash: transient reactance of d-axis, given in [pu]
- X_q_dash: transient reactance of q-axis, given in [pu]
- X_d: reactance of d-, given in [pu]
- X_d: reactance of q-axis, given in [pu]

IEEE DC1A Exciter Model
- K_e: Exciter constant related to self-excited field [pu]
- K_f: Excitation control system stabilizer gains [pu]
- K_a: Voltage Regulator gain [pu]
- U:
- U_ref: Reference value of the stator terminal voltage [pu]
- U_ref2: Reference value of the stator terminal voltage [pu]
- U_rmax: Voltage regulator maximum output [pu]
- U_rmin: Voltage regulator minimum output [pu]
- T_a: Time constant of the voltage regulator [s]
- T_f: Excitation control system stabilizer time constant [s]
- T_e: Exciter time constant, integration rate associated with exciter control [s]

"""

using PowerDynamics
using NetworkDynamics
import PowerDynamics: dimension, symbolsof, construct_vertex
import PowerDynamics: AbstractNode
import Base: @__doc__
include("ExciterSaturtionEq.jl")

@DynamicNode ThirdOrderEqExciterIEEED1A(H, P, D, Ω, T_d_dash, X_d_dash, X_q_dash, X_d, K_e, K_f, K_a, U, U_ref, U_ref2, U_rmax, U_rmin, T_a, T_f, T_e, S_E_max, S_E_tq, V_R_max) begin
    @assert H > 0 "inertia (H) should be >0"
    @assert D >= 0 "damping (D) should be >=0"
    @assert T_d_dash > 0 "time constant of d-axis (T_d_dash) should be >0"
    @assert X_d_dash >= 0 "transient reactance of d-axis (X_d_dash) should be >=0"
    @assert X_d >= 0 "reactance of d-axis (X_d_dash) should be >=0"

    Ω_H = Ω / (2 * H)
end [[θ,dθ],[ω, dω], [E_f, dE_f], [U_f, dU_f], [U_r, dU_r]] begin
    i_c = 1im * i * exp(-1im * θ)
    e_c = 1im * u * exp(-1im * θ)
    p = real(u * conj(i))
    e_d = 0           # simpilfication from K.Schmietendorf paper
    e_q = imag(e_c)
    i_d = real(i_c)
    i_q = imag(i_c)

    # Exciter Model, IEEEDC1A
	E_f_d_max, A_x, B_x = ExciterSaturationEq(K_e, S_E_max, S_E_tq, V_R_max) # Saturation modeling
	U_x = A_x * exp(B_x * E_f)

	dU_r = 1 / T_a * (K_a * (U_ref - U + U_ref2 - U_f) - U_r)
	dU_f = 1 / T_f * (K_f / T_e * (U_r - U_x - K_e * E_f) - U_f)

    if U_r > U_rmax
		println("Ur > Urmax")
		U_r2 = U_rmax

	elseif U_r < U_rmin
		println("Ur < Urmax")
		U_r2 = U_rmin
	else
		U_r2 = U_r
	end

	dE_f = 1 / T_e * ( U_r2 - U_x - K_e * E_f)

    # Third Order Model
	dθ = ω
    dω = (P - D * ω - p - (X_q_dash - X_d_dash) * i_d * i_q) * Ω_H
    de_q = (1 / T_d_dash) * (E_f - e_q + i_d * (X_d - X_d_dash))
    # -> u = e_q * exp(1im * θ)
    du = de_q * exp(1im * θ) + u * 1im * ω

end

export ThirdOrderEqExciterIEEEDoneA
