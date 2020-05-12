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
append!(line_list, [StaticLine(from = 1, to = 3, Y = -1im / 0.02)])


node_list = []
append!(node_list, [SlackAlgebraic(U = 1, Y_n = 1)])
append!(node_list, [ThirdOrderEq(H = 5.0, P = 1.0, D = 0.1, Ω = 50, E_f = 0.2, X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1)])
append!(node_list, [ThirdOrderEq(H = 5.0, P = 1.0, D = 0.1, Ω = 50, E_f = 0.2,X_d_dash = 0.111, T_d_dash = 0.1, X_q_dash = 0.103, X_d = 0.1)])


powergrid = PowerGrid(node_list, line_list)

operationpoint = find_operationpoint(powergrid)

result = simulate(
    Perturbation(2, :ω, Inc(0.1)),
    powergrid,
    operationpoint,
    timespan = (0.0, 10),
)

plot_res(result, powergrid, 2)
