using PowerDynamics
using NetworkDynamics
using LaTeXStrings
include("plotting.jl")
include("ThirdOrderEq.jl")
include("ThirdOrderEqGovernorIEEEG1.jl")
include("ThirdOrderEqExciterIEEEDC1A.jl")

line_list = []
append!(line_list, [StaticLine(from = 1, to = 2, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 3, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 3, to = 4, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 2, to = 4, Y = -1im / 0.02)])
append!(line_list, [StaticLine(from = 1, to = 4, Y = -1im / 0.02)])



node_list = []
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 5.0, P = 0.6, D = 0.01, Ω = 50, E_f = 0.2, X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1)])
append!(node_list, [ThirdOrderEq(H = 5.0, P = -0.6, D = 0.01, Ω = 50, E_f = 0.2, X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1)])
append!(node_list, [ThirdOrderEqExciterIEEED1A(S_E_max = 0.86, S_E_tq = 0.5, V_R_max = 7.6,H = 5.0, P = 1.0, D = 0.01,  Ω = 50,U = 0.2, X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1, K_e = 0.2 , K_f = 0.4, K_a = 0.5, U_ref = 1.2, U_ref2 = 0.9, U_rmax = 1.75, U_rmin = 0.0, T_e = 0.8, T_f = 0.3, T_a = 0.3)])

powergrid = PowerGrid(node_list, line_list)

operationpoint = find_operationpoint(powergrid)

result = simulate(
    Perturbation(3, :ω, Inc(0.00001)),
    powergrid,
    operationpoint,
    timespan = (0.0, 50),
)

plot_res(result, powergrid, 3)
