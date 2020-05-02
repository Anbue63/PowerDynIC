using Revise
using PowerDynamics
using Plots
using DifferentialEquations
import PowerDynamics:rhs
# some additional functions that might develop into new PD.jl features
include("PDpatches.jl")

# custom node type
include("components/ThirdOrderEq.jl")

# load example system
include("../example_cases/BasinStabilityForPD/system.jl")

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end

powergrid = PowerGrid(node_list,line_list)
rpg = rhs(powergrid)

state(vec) = State(powergrid, vec)
project(vec) = find_valid_initial_condition(powergrid, vec)


## 1 IC guess

ic1 = ones(systemsize(powergrid)) |> state
ic0 = zeros(systemsize(powergrid)) |> state
icrand = rand(systemsize(powergrid)) |> state

# let's try a type-specific guess
sl = node_list[1]
guess(::Type{SlackAlgebraic}) = [sl.U, 0.] #[:u_r, :u_i]
guess(::Type{PQAlgebraic}) = [sl.U, 0.] #[:u_r, :u_i]
guess(::Type{ThirdOrderEq}) = [sl.U, 0., 0., 0.] #[:u_r, :u_i, :θ, :ω]
type_guesses = node_list .|> typeof .|> guess
icm = vcat(type_guesses...) |> state

op_ic1 = find_operationpoint(powergrid; ic_guess=ic1.vec, tol=1E-10)
op_ic0 = find_operationpoint(powergrid; ic_guess=ic0.vec, tol=1E-10)
op_icm = find_operationpoint(powergrid; ic_guess=icm.vec, tol=1E-10)

# random is random ...
op_icrand = find_operationpoint(powergrid; ic_guess=icrand.vec, tol=1E-14)


# we can try to improve nlsolve be using the sparse structure
# this only works for icm, the others fail
sp_op_icm = find_operationpoint_sparse(powergrid; ic_guess=icm.vec, tol=1E-10)

begin
    scatter(op_ic1[:, :u], label="1", legend=:topleft)
    scatter!(op_ic0[:, :u], label="0")
    scatter!(op_icrand[:, :u], label="rand")
    scatter!(op_icm[:, :u], label="mix")
    scatter!(sp_op_icm[:, :u], label="mix sp")
end

## 2 use ForwardDiff to check stability

# apparently the op's are not stable
# actually I am not sure whether eigvals(j(vec) * pinv(M) * M) is correct
check_eigenvalues(rpg, op_ic1)
check_eigenvalues(rpg, op_ic0)
check_eigenvalues(rpg, op_icrand)

# positive EV but small, maybe actually stable
check_eigenvalues(rpg, op_icm)

# this one seems to be stable!
check_eigenvalues(rpg, sp_op_icm)

## 3 compare with SteadyStateProblem in DiffEq

perturb = Perturbation(7, :ω, Inc(0.1))

function find_steady_state(rpg, u0)
    ode = ODEProblem(rpg, u0.vec, (0., 20.), nothing)
    ssp = SteadyStateProblem(ode)
    return solve(ssp, DynamicSS(Rodas5(); abstol=1e-8,reltol=1e-6,tspan=Inf)) |> Array |> state
end

de_sp_op_icm = find_steady_state(rpg, perturb(sp_op_icm))

begin
    scatter(sp_op_icm[:, :u], ms=10, label="mix sp", legend=:topleft)
    scatter!(de_sp_op_icm[:, :u], label="SSP")
end

# is this method maybe more robust?
de_op_icm = find_steady_state(rpg, perturb(op_icm))
de_op_ic1 = find_steady_state(rpg, perturb(op_ic1))
# this one diverges
#de_op_ic0 = find_steady_state(rpg, perturb(op_ic0))
de_op_icrand = find_steady_state(rpg, perturb(op_icrand))

begin
    scatter(de_op_ic1[:, :u], label="1", legend=:topleft)
    scatter!(de_op_icrand[:, :u], label="rand")
    scatter!(de_op_icm[:, :u], label="mix")
    scatter!(de_sp_op_icm[:, :u], label="mix sp")
end

# they are all basically stable within 1E-6
check_eigenvalues(rpg, de_op_ic1)
check_eigenvalues(rpg, de_op_icm)
check_eigenvalues(rpg, de_op_icrand)
check_eigenvalues(rpg, de_sp_op_icm)

## 4 check integration as ODE problem with mass_matrix

ω_idx = findall(:ω .∈ symbolsof.(node_list))

function sim(rpg, u0)
    ode = ODEProblem(rpg, u0.vec, (0., 20.), nothing)
    dqsol = solve(ode, Rodas5(), reltol=1e-12, abstol=1e-12)
    return PowerGridSolution(dqsol, powergrid)
end

# -> looks stable :-)
perturb = Perturbation(7, :ω, Inc(0.1))
sol = sim(rpg, perturb(de_sp_op_icm))
plot(sol, ω_idx, :ω, legend=false)


# now try the case in test_basinstability
# -> it gives dtmin error and stops
perturb = Perturbation(7, :ω, Inc(29.999995))

sol = sim(rpg, perturb(de_sp_op_icm))
plot(sol, ω_idx, :ω, legend=false)

## 6 check Sundials solvers
# https://docs.sciml.ai/stable/solvers/dae_solve/#Sundials.jl-1

using Sundials

ω_idx = findall(:ω .∈ symbolsof.(node_list))

#rpg_dense = ODEFunction(rpg.f, mass_matrix=Array(rpg.mass_matrix))

function sim(rpg, u0)
    ode = ODEProblem(rpg, u0.vec, (0., 20.), nothing)
    dqsol = solve(ode, IDA(linear_solver=:GMRES), reltol=1e-12, abstol=1e-12)
    return PowerGridSolution(dqsol, powergrid)
end

perturb = Perturbation(7, :ω, Inc(0.1))
sol = sim(rpg, perturb(de_sp_op_icm))
plot(sol, ω_idx, :ω, legend=false)

## 7 investigate DAE formulation

using LinearAlgebra
using Sundials

# load example system
include("../example_cases/DAE_Solver_PD/system.jl")

powergrid = PowerGrid(node_list,line_list)

state(vec) = State(powergrid, vec)
ω_idx = findall(:ω .∈ symbolsof.(node_list))

# ND.jl gives us an ODEFunction so we have to construct one manually
rpg_ode = rhs(powergrid)

function dae_form(res, du, u, p, t)
    du_temp = similar(du)
    rpg_ode.f(du_temp, u, p, t)
    @. res = du - du_temp
end

# construct initial guess
sl = node_list[1]
guess(::Type{SlackAlgebraic}) = [sl.U, 0.] #[:u_r, :u_i]
guess(::Type{PQAlgebraic}) = [sl.U, 0.] #[:u_r, :u_i]
guess(::Type{ThirdOrderEq}) = [sl.U, 0., 0., 0.] #[:u_r, :u_i, :θ, :ω]
type_guesses = node_list .|> typeof .|> guess
icm = vcat(type_guesses...) |> state

# find PD operation point
op_icm = find_operationpoint(powergrid; ic_guess=icm.vec, tol=1E-10)

u0 = op_icm.vec
du0 = similar(u0)
diff_vars = (diag(rpg_ode.mass_matrix) .== 1)
p = nothing
tspan = (0., 20.)

rpg_dae = DAEFunction{true, true}(dae_form, syms=rpg_ode.syms)
dae = DAEProblem(rpg_dae, du0, u0, tspan, p; differential_vars=diff_vars)

#ToDO : test other solvers, e.g. IDA(linear_solver=:GMRES)
# ShampineCollocationInit
sol = solve(dae, IDA(linear_solver=:GMRES), initializealg=BrownFullBasicInit())

# apparently, this doesn't run into a stable state
plot(PowerGridSolution(sol, powergrid), ω_idx, :ω, legend=false)
plot(PowerGridSolution(sol, powergrid), :, :v, legend=false)

perturb = Perturbation(2, :ω, Inc(31.))
u0 = perturb(op_icm)
du0 = similar(u0.vec)

dae = DAEProblem(rpg_dae, du0, u0.vec, (0., 100.); differential_vars=diff_vars)
dqsol = solve(dae, IDA(linear_solver=:GMRES), initializealg=BrownFullBasicInit())
pgsol = PowerGridSolution(dqsol, powergrid)

plot(pgsol, ω_idx, :ω, legend=false)
plot(pgsol, :, :v, legend=false)
