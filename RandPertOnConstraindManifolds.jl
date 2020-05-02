using PowerDynamics
using ForwardDiff
using ForwardDiff: Dual
using Distributions:Uniform
import PowerDynamics:rhs
using DifferentialEquations
using LinearAlgebra
using DiffEqCallbacks
include("system.jl")
include("/home/anna/Documents/MASTER_ARBEIT/PowerDynIC-master/src/PDpatches.jl")
include("/home/anna/Documents/MASTER_ARBEIT/PowerDynIC-master/src/PDpatches.jl")

"""
A method to find random pertubations in the case of systems with algebraic
constraints. There is then a set of valid initial conditions which lies on a so
called constrained manifold.
The system must still lie on this manifold after the random pertubation.
This could be applied to Basin Stability Methods in the case of complexer
Generator Methods which contain algebraic constraints.
"""
function RandPertWithConstrains(pg::PowerGrid, operationpoint,node, interval)
    g(θ,i,u) = [1im * i * exp(-1im * θ), 1im * u * exp(-1im * θ)] # Constraints

    z = [operationpoint[node,:θ],operationpoint[node,:iabs],operationpoint[node,:v]]

    Jacg_imag = ForwardDiff.jacobian(x -> imag(g(x[1],x[2],x[3])), z)
    Jacg_real = ForwardDiff.jacobian(x -> real(g(x[1],x[2],x[3])), z)

    Jacg = Jacg_real + 1im .* Jacg_imag  # I could not find a package which is able to handle complex valued function

    kerJacg = nullspace(Jacg) # kernel of Jacg, returns the extrem points/saddle points

    Pro_kerg = Projvw(z, kerJacg) # Projection of z on ker(∇g(z))

    Frand = rand(Uniform(interval[1],interval[2])) # random pertubation "force"

    dz = Pro_kerg * Frand # Random Pertubation along the Projection of the Constraind Manifold

    PerturbedState = copy(operationpoint) # otherwise operationpoint is overwritten
    PerturbedState[node, :θ] = abs(dz[1]) # thats not correct but its a quick workaround for now
    PerturbedState[node, :v] = dz[3]

    return PerturbedState
end

#=
PerturbedState = RandPertWithConstrains(powergrid, op, 3, [-100,100])

result = simulate(Perturbation(3,:v,Inc(0)), powergrid, PerturbedState,timespan = 20.0)

ω_idx = findall(:ω .∈ symbolsof.(powergrid.nodes))

θ_idx = findall(:θ .∈ symbolsof.(powergrid.nodes))
plot(PowerGridSolution(result.dqsol, powergrid), θ_idx, :θ, legend=false)
plot(PowerGridSolution(result.dqsol, powergrid), ω_idx, :ω, legend=false)

plot(PowerGridSolution(result.dqsol, powergrid), :, :v, legend=false)

=#
"""
Calculates and returns the Projection of a vector w on the vector v.
"""
function Projvw(v,w)
    dot(w,v) * v * (1 / norm(v)^2)
end

function rhs(pg::PowerGrid)
    network_dynamics(map(construct_vertex, pg.nodes), map(construct_edge, pg.lines), pg.graph)
end



"""
A method to assures that the conditions given from g stay on the contrained
manifold. A DifferentialEquations Callback is used.
This method is not working yet. TODO: I need to find out how PD fetches the equations
from the Node Makro. I will then use this algorithm to feed in the eqs here.
"""
function SimManifoldCd(pg::PowerGrid, u0::State, node,timespan)

    ode = ODEProblem(rhs(pg), u0.vec, (0., timespan), nothing)

    function g(resid,sym,p,t)
        resid[1] = 1im * sym[2] * exp(-1im * sym[3])
        resid[2] = 1im * sym[3] * exp(-1im * sym[3])
    end
    cb = ManifoldProjection(g)


    println("RIGHTHANDSIDE:    ",rhs(pg))


    sol = solve(ode,Rodas5(),save_everystep=false,callback=cb, abstol=1e-12)
    dqsol = solve(ode, Rodas5(), reltol=1e-12, abstol=1e-12)
    return PowerGridSolution(dqsol, powergrid)
end
