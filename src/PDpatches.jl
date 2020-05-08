import PowerDynamics: rhs, get_current
using DifferentialEquations: DAEFunction
using NetworkDynamics: network_dynamics, StaticEdgeFunction

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph; parallel=true)
end

get_current(state, n) = begin
    vertices = map(construct_vertex, state.grid.nodes)
    edges = map(construct_edge, state.grid.lines)
    sef = StaticEdgeFunction(vertices,edges,state.grid.graph, parallel=true)
    (e_s, e_d) = sef(state.vec, Nothing, 0)
    total_current(e_s[n], e_d[n])
end

function rhs_dae(pg::PowerGrid)
    rpg_ode = rhs(pg)

    du_temp = zeros(systemsize(pg))

    function dae_form(res, du, u, p, t)
        rpg_ode.f(du_temp, u, p, t)
        @. res = du - du_temp
    end

    return DAEFunction{true, true}(dae_form, syms=rpg_ode.syms)
end

using LinearAlgebra: diag
using OrdinaryDiffEq: ODEFunction

function differential_vars(pg::PowerGrid)
    rpg_ode = rhs(pg)
    return (diag(Array(rpg_ode.mass_matrix)) .== 1)
end

function differential_vars(rpg_ode::ODEFunction)
    return (diag(Array(rpg_ode.mass_matrix)) .== 1)
end

import PowerDynamics: find_operationpoint, State
using NLsolve: nlsolve, converged

function find_operationpoint(pg::PowerGrid; ic_guess = nothing, tol=1E-9)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find the fixed point"
    end
    if SwingEq ∈ pg.nodes .|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end

    if ic_guess === nothing
        ic_guess = initial_guess(pg)
    end

    rr = RootRhs(rhs(pg))
    nl_res = nlsolve(rr, ic_guess; xtol = tol)
    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        throw(OperationPointError("Failed to find initial conditions on the constraint manifold!"))
    end
end

using ForwardDiff: jacobian
using SparseArrays: sparse
using LinearAlgebra: eigvals
using NLsolve: OnceDifferentiable

function check_eigenvalues(rpg, vec)
    M = Array(rpg.mass_matrix)
    f!(dx, x) = rpg(dx, x, nothing, 0.)
    j(x) = (dx = similar(x); jacobian(f!, dx, x))
    λ = eigvals(j(vec) * pinv(M) * M) .|> real |> extrema
    println("Jacobian spectrum \nmin : ", first(λ), "\nmax : ",last(λ), "\nstable : ", isapprox(last(λ), 0, atol=1e-8))
    return λ
end

check_eigenvalues(rpg, s::State) = check_eigenvalues(rpg, s.vec)

function find_operationpoint_sparse(pg::PowerGrid; ic_guess = nothing, tol=1E-9)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers. Currently not making any checks concerning assumptions of whether its possible to find the fixed point"
    end
    if SwingEq ∈ pg.nodes .|> typeof
        throw(OperationPointError("found SwingEq node but these should be SwingEqLVS (just SwingEq is not yet supported for operation point search)"))
    end

    if ic_guess === nothing
        system_size = systemsize(pg)
        ic_guess = ones(system_size)
    end

    # we can specify the Jacobian to be sparse
    rpg = rhs(pg)
    f!(dx, x) =  rpg(dx, x, nothing, 0.)
    f(x) = (dx = similar(x); f!(dx, x); dx) # this is the same as RootRhs
    j!(J, x) = (J .= jacobian(f, x))

    n = systemsize(pg)

    dx0 = similar(ic_guess)
    F0 = similar(ic_guess)
    J0 = sparse(zeros(n, n))
    df = OnceDifferentiable(f!, j!, dx0, F0, J0)
    nl_res = nlsolve(df, ic_guess; xtol = tol)

    if converged(nl_res) == true
        return State(pg, nl_res.zero)
    else
        throw(OperationPointError("Failed to find initial conditions on the constraint manifold!"))
    end
end

using OrdinaryDiffEq: ODEFunction
using NLsolve: nlsolve, converged
using LinearAlgebra: pinv

struct RootRhs_ic
    rhs
    mpm
end
function (rr::RootRhs_ic)(x)
    dx = similar(x)
    rr.rhs(dx, x, nothing, 0.)
    rr.mpm * dx .- dx
end

function RootRhs_ic(of::ODEFunction)
    mm = of.mass_matrix
    @assert mm != nothing
    mpm = pinv(mm) * mm
    RootRhs_ic(of.f, mpm)
end


function find_valid_initial_condition(pg::PowerGrid, ic_guess)
    rr = RootRhs_ic(rhs(pg))
    nl_res = nlsolve(rr, ic_guess)
    if converged(nl_res) == true
        return nl_res.zero #State(pg, nl_res.zero)
    else
        println("Failed to find initial conditions on the constraint manifold!")
        println("Try running nlsolve with other options.")
    end
end

using PowerDynamics: AbstractNode
"""
Makes an type specific initial guess to help the operation point search.
The voltage is of all nodes is fixed to the voltage of the first SlackAlgebraic
in the system. The other symbols are set to zero.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    guess: Type specific initial guess
"""
function initial_guess(pg)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers."
    end

    sl = findfirst(SlackAlgebraic  ∈  pg.nodes .|> typeof)
    slack = pg.nodes[sl]

    type_guesses = guess.(pg.nodes, slack.U)
    return vcat(type_guesses...)
end


"""
guess(::Type{SlackAlgebraic}) = [slack.U, 0.]         #[:u_r, :u_i]
guess(::Type{PQAlgebraic}) = [slack.U, 0.]            #[:u_r, :u_i]
guess(::Type{ThirdOrderEq}) = [slack.U, 0., 0., 0.]   #[:u_r, :u_i, :θ, :ω]
"""
function guess(n::AbstractNode, slack_voltage)
    voltage = zeros(dimension(n))
    voltage[1] = real(slack_voltage)  #[:u_r, :u_i]
    voltage[2] = imag(slack_voltage)
    return voltage
end

function guess(n::SlackAlgebraic, slack_voltage)
    @assert n.U == slack_voltage
    return [real(n.U), imag(n.U)]
end

using OrdinaryDiffEq: ODEFunction
using LightGraphs: AbstractGraph
# stack dynamics with higher level controller
struct ControlledPowerGrid
    graph:: G where G <: AbstractGraph
    nodes
    lines
    controller # oop function
end

# Constructor
function ControlledPowerGrid(control, pg::PowerGrid, args...)
    ControlledPowerGrid(pg.graph, pg.nodes, pg.lines, (u, p, t) -> control(u, p, t, args...))
end

function rhs(cpg::ControlledPowerGrid)
    open_loop_dyn = PowerGrid(cpg.graph, cpg.nodes, cpg.lines) |> rhs
    function cpg_rhs!(du, u, p, t)
        # Get controller state
        p_cont = cpg.controller(u, p, t)
        # Calculate the derivatives, passing the controller state through
        open_loop_dyn(du, u, p_cont, t)
    end
    ODEFunction{true}(cpg_rhs!, mass_matrix=open_loop_dyn.mass_matrix, syms=open_loop_dyn.syms)
end
