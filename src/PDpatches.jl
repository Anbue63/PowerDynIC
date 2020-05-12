import PowerDynamics: rhs
#=
function rhs(pg::PowerGrid; parallel=false)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph; parallel=parallel)
end
=#
import PowerDynamics: find_operationpoint
using NLsolve: nlsolve, converged

function find_operationpoint(pg::PowerGrid; ic_guess = nothing, tol=1E-9)
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

"""
Ist das so richtig beschrieben?
Calculates and returns the eigenvalues λ of the jacobian of the right hand side of
the powergrid at the state s.
The sign of the real part of the eigenvalues will decide whether a fixed point
will be attracting (λ_max < 0) or repelling (λ_max > 0).
Inputs:
       pg: Power grid, a graph containg nodes and lines
       s: State of the powergrid
Outputs:
        λ: The biggest and the smallest eigenvalue
        stable: Determines the stability by checking if λ_max is positive
"""
function check_eigenvalues(pg::PowerGrid, s::State)
    rpg = rhs(pg)
    M = Array(rpg.mass_matrix)
    f!(dx, x) = rpg(dx, x, nothing, 0.)
    j(x) = (dx = similar(x); jacobian(f!, dx, x))
    λ = eigvals(j(s.vec) * pinv(M) * M) .|> real |> extrema
    stable = isapprox(last(λ), 0, atol=1e-8)
    return λ, stable
end

#check_eigenvalues(rpg, s::State) = check_eigenvalues(rpg, s.vec)

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

"""
Makes an type specific initial guess to help the operation point search.
The voltage is of all nodes is fixed to the voltage of the first SlackAlgebraic
in the system. The other symbols are set to zero.
Inputs:
    pg: Power grid, a graph containg nodes and lines
Outputs:
    guess: Type specific initial guess
"""
function initial_guess(pg::PowerGrid)
    if SlackAlgebraic ∉ pg.nodes .|> typeof
        @warn "There is no slack bus in the system to balance powers."
    end

    sl = findfirst(SlackAlgebraic  ∈  pg.nodes .|> typeof)
    slack = pg.nodes[sl]
    guess(::Type{SlackAlgebraic}) = [slack.U, 0.]         #[:u_r, :u_i]
    guess(::Type{PQAlgebraic}) = [slack.U, 0.]            #[:u_r, :u_i]
    guess(::Type{ThirdOrderEq}) = [slack.U, 0., 0., 0.]   #[:u_r, :u_i, :θ, :ω]

    type_guesses = pg.nodes .|> typeof .|> guess
    icm = vcat(type_guesses...)
end
